"""Analyze an audited promoter/CDS bundle; never infer missing sample identities."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import sys

import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle, Wedge
import networkx as nx
import numpy as np
import pandas as pd
from scipy import stats

from common import MISSING, empty_output, reconstruct_proteins, sha256, write_json
from style import (BASE_COLORS, HAP_COLORS, POPULATIONS, POP_COLORS, POP_LABELS,
                   RESPONSE_COLORS, save, setup)


def read(path, columns):
    data = pd.read_csv(path, dtype=str, keep_default_na=False)
    if not set(columns) <= set(data):
        raise ValueError(f'{path.name}: required columns {columns}')
    return data


def validate(haps, sites, counts, assignments=None):
    ids = haps.haplotype_id.tolist()
    keys = sites.site_id.tolist()
    if len(ids) != len(set(ids)) or len(keys) != len(set(keys)):
        raise ValueError('Duplicate haplotype or site identifiers')
    if not set(keys) <= set(haps):
        raise ValueError('Haplotype matrix lacks declared sites')
    if keys and haps.duplicated(keys).any():
        raise ValueError('Duplicate allele combinations must be merged before analysis')
    if 'source_id' in haps and haps.source_id.duplicated().any():
        raise ValueError('Duplicate source haplotype identifiers')
    positions = pd.to_numeric(sites.position, errors='raise')
    if (positions < 1).any() or (positions % 1 != 0).any():
        raise ValueError('Genomic positions must be positive 1-based integers')
    for row in sites.to_dict('records'):
        allowed = str(row['alleles']).upper().split('/')
        for value in haps[row['site_id']]:
            if str(value).upper() not in allowed and str(value).upper() not in MISSING:
                raise ValueError(f'Undeclared allele at site {row["site_id"]}: {value}')
    if counts.duplicated(['haplotype_id', 'population']).any():
        raise ValueError('Duplicate count rows')
    expected = {(h, p) for h in ids for p in POPULATIONS}
    if set(zip(counts.haplotype_id, counts.population)) != expected:
        raise ValueError('Counts must contain each haplotype × population exactly once')
    numbers = pd.to_numeric(counts['count'], errors='raise')
    if not np.isfinite(numbers).all() or (numbers < 0).any() or (numbers % 1 != 0).any():
        raise ValueError('Counts must be nonnegative integers')
    counts['count'] = numbers.astype(int)
    if assignments is not None:
        if not set(assignments.haplotype_id) <= set(ids) or not set(assignments.population) <= set(POPULATIONS):
            raise ValueError('Unknown sample assignment')
        groups = assignments.groupby('sample_id')
        if (groups.size() > 2).any() or (groups.population.nunique() > 1).any():
            raise ValueError('Assignments allow at most two copies and one population per material')
        observed = assignments.groupby(['haplotype_id', 'population']).size()
        for row in counts.itertuples():
            if observed.get((row.haplotype_id, row.population), 0) != row.count:
                raise ValueError('Sample copy assignments do not reproduce summary counts')


def materials(assignments):
    """Collapse homozygous copies; exclude ambiguous/heterozygous material."""
    if assignments is None:
        return None, pd.DataFrame(columns=['sample_id', 'reason'])
    ntypes = assignments.groupby('sample_id').haplotype_id.nunique()
    bad = ntypes[ntypes > 1].index
    excluded = pd.DataFrame({'sample_id': bad, 'reason': 'multiple retained haplotypes'})
    unique = assignments[~assignments.sample_id.isin(bad)].drop_duplicates('sample_id')
    return unique, excluded


def association(unique, label):
    result = {'level': label, 'unit': 'independent material', 'p': None}
    if unique is None or unique.empty:
        return {**result, 'status': 'skipped: no unique material assignments'}
    table = pd.crosstab(unique.haplotype_id, unique.population)
    table = table.loc[table.sum(axis=1) > 0, table.sum(axis=0) > 0]
    result['n'] = int(table.to_numpy().sum())
    if min(table.shape) < 2:
        return {**result, 'status': 'skipped: fewer than two observed groups'}
    obs = table.to_numpy()
    chi, p, _, expected = stats.chi2_contingency(obs, correction=False)
    result.update({'statistic_chi2': float(chi),
                   'cramers_v': float(np.sqrt(chi / (obs.sum() * (min(obs.shape) - 1))))})
    if (expected >= 5).all():
        method = 'Pearson chi-square'
    elif obs.shape == (2, 2):
        _, p = stats.fisher_exact(obs)
        method = 'Fisher exact, two-sided'
    else:
        rng = np.random.default_rng(2026)
        simulated = stats.random_table.rvs(obs.sum(axis=1), obs.sum(axis=0), size=9999, random_state=rng)
        score = ((simulated - expected) ** 2 / expected).sum(axis=(1, 2))
        p = (np.count_nonzero(score >= chi - 1e-12) + 1) / 10000
        method = 'fixed-margin Monte Carlo chi-square; 9999 draws; seed 2026'
    return {**result, 'status': 'ok', 'method': method, 'p': float(p)}


def adjust_bh(results):
    valid = sorted([r for r in results if r['p'] is not None], key=lambda r: r['p'])
    running = 1.0
    for i in range(len(valid) - 1, -1, -1):
        running = min(running, valid[i]['p'] * len(valid) / (i + 1))
        valid[i]['q_bh_this_run'] = running


def matrix_plot(haps, keys, labels, out, name, title, amino=False):
    if not keys or haps.empty:
        return
    widths = [max(1, max(len(str(v)) for v in haps[k]) * .32 + .5) for k in keys]
    edges = np.r_[0, np.cumsum(widths)]
    centers = (edges[:-1] + edges[1:]) / 2
    fig, ax = plt.subplots(figsize=(max(6, sum(widths) * .65 + 2), max(3, len(haps) * .5 + 2)))
    states = sorted({str(x) for k in keys for x in haps[k]})
    palette = {v: HAP_COLORS[i % len(HAP_COLORS)] for i, v in enumerate(states)}
    for y, row in enumerate(haps.to_dict('records')):
        for x, key in enumerate(keys):
            value = str(row[key])
            color = palette[value] if amino else '#EEEEEE' if value.upper() in MISSING else BASE_COLORS.get(value, '#BBBBBB')
            ax.add_patch(Rectangle((edges[x], y - .5), widths[x], 1, facecolor=color, edgecolor='white'))
            ax.text(centers[x], y, value, ha='center', va='center', fontsize=12,
                    color='white' if value == 'T' and not amino else 'black')
    ax.set(xlim=(0, edges[-1]), ylim=(len(haps) - .5, -.5), title=title)
    ax.set_xticks(centers, labels, rotation=60, ha='right')
    ax.set_yticks(range(len(haps)), haps.haplotype_id)
    ax.set_xlabel('Amino acid position (1-based)' if amino else 'Retained variant position')
    ax.set_ylabel('Haplotype')
    save(fig, out, name)


def frequencies(counts, ids, out, prefix):
    table = counts.pivot(index='haplotype_id', columns='population', values='count').reindex(index=ids, columns=POPULATIONS).fillna(0)
    pop_total, hap_total = table.sum(axis=0), table.sum(axis=1)
    source = counts.copy()
    source['population_total_copies'] = source.population.map(pop_total)
    source['haplotype_total_copies'] = source.haplotype_id.map(hap_total)
    source['pct_within_population'] = source['count'].div(source.population_total_copies.replace(0, np.nan)) * 100
    source['pct_within_haplotype'] = source['count'].div(source.haplotype_total_copies.replace(0, np.nan)) * 100
    source.to_csv(out / f'{prefix}_frequencies.csv', index=False)
    for within_pop in [True, False]:
        labels = POPULATIONS if within_pop else ids
        series = ids if within_pop else POPULATIONS
        vals = table.T if within_pop else table
        denominator = pop_total if within_pop else hap_total
        fig, ax = plt.subplots(figsize=(max(7, len(labels) * .7 + 3), 4.8))
        bottom = np.zeros(len(labels))
        for i, category in enumerate(series):
            values = vals[category].to_numpy() / denominator.replace(0, np.nan).to_numpy() * 100
            ax.bar(range(len(labels)), values, bottom=bottom,
                   color=HAP_COLORS[i % len(HAP_COLORS)] if within_pop else POP_COLORS[category],
                   label=category if within_pop else POP_LABELS[category])
            bottom += np.nan_to_num(values)
        for i, n in enumerate(denominator):
            ax.text(i, 102, f'N={n}' if n else 'No observations', ha='center', fontsize=12)
        ax.set_xticks(range(len(labels)), [POP_LABELS.get(x, x) for x in labels], rotation=0 if within_pop else 45)
        ax.set(ylim=(0, 115), ylabel='Percent of returned allele copies (%)',
               title='Haplotypes within each population' if within_pop else 'Populations within each haplotype')
        ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', frameon=False)
        save(fig, out / 'figures', f'{prefix}_within_population' if within_pop else f'{prefix}_within_haplotype')
    return table


def network(haps, keys, table, out):
    if not keys or any(str(v).upper() in MISSING for k in keys for v in haps[k]):
        return 'skipped: missing states or no retained sites'
    ids = haps.haplotype_id.tolist()
    matrix = haps[keys].to_numpy()
    graph = nx.Graph()
    graph.add_nodes_from(ids)
    for i, a in enumerate(ids):
        for j in range(i):
            graph.add_edge(a, ids[j], weight=int(np.count_nonzero(matrix[i] != matrix[j])))
    tree = nx.minimum_spanning_tree(graph)
    pd.DataFrame([{'source': a, 'target': b, 'site_differences': d['weight']} for a, b, d in tree.edges(data=True)],
                 columns=['source', 'target', 'site_differences']).to_csv(out / 'network_edges.csv', index=False)
    table.reset_index().to_csv(out / 'network_nodes.csv', index=False)
    # A layered tree drawing avoids crossing edges. The anchor is not an ancestor.
    anchor = nx.center(tree)[0]
    pos, cursor = {}, [0]

    def place(node, parent=None, depth=0):
        children = [n for n in tree.neighbors(node) if n != parent]
        child_x = [place(n, node, depth + 1) for n in children]
        x = float(np.mean(child_x)) if child_x else cursor[0] * 2.5
        if not child_x:
            cursor[0] += 1
        pos[node] = (x, -depth * 2)
        return x

    place(anchor)
    depth = max(-p[1] / 2 for p in pos.values())
    fig, ax = plt.subplots(figsize=(max(8, cursor[0] * 2 + 3), max(6, depth * 1.3 + 2)))
    nx.draw_networkx_edges(tree, pos, ax=ax, edge_color='#777777')
    nx.draw_networkx_edge_labels(tree, pos, edge_labels=nx.get_edge_attributes(tree, 'weight'), font_size=12, rotate=False, ax=ax)
    max_n = max(float(table.sum(axis=1).max()), 1)
    for hid, (x, y) in pos.items():
        values = table.loc[hid]
        radius = .35 * np.sqrt(values.sum() / max_n)
        if radius == 0:
            ax.scatter([x], [y], marker='x', color='gray')
        angle = 0
        for pop in POPULATIONS:
            theta = 360 * values[pop] / values.sum() if values.sum() else 0
            if theta:
                ax.add_patch(Wedge((x, y), radius, angle, angle + theta, facecolor=POP_COLORS[pop], edgecolor='white', linewidth=.5))
            angle += theta
        ax.text(x + .55, y, f'{hid}\nN={int(values.sum())}', ha='left', va='center', fontsize=12)
    ax.legend(handles=[Patch(color=POP_COLORS[p], label=POP_LABELS[p]) for p in POPULATIONS],
              bbox_to_anchor=(1, 1), loc='upper left', frameon=False)
    ax.set(xlim=(-1, max(p[0] for p in pos.values()) + 2.5), ylim=(-depth * 2 - 1, 1), aspect='equal',
           title='Retained-site minimum spanning tree\nEdge labels: differing site states; node area: allele copies')
    ax.axis('off')
    save(fig, out / 'figures', 'haplotype_network_mst')
    return 'ok: minimum spanning tree; no ancestry inferred'


def geography(unique, metadata_path, out):
    if unique is None or not metadata_path.exists():
        return 'skipped: material assignments or coordinates unavailable'
    meta = read(metadata_path, ['sample_id', 'longitude', 'latitude'])
    if meta.sample_id.duplicated().any():
        raise ValueError('Duplicate material metadata; reconcile before joining')
    data = unique.merge(meta, on='sample_id', how='left', validate='one_to_one', suffixes=('', '_metadata'))
    for key in ['longitude', 'latitude']:
        data[key] = pd.to_numeric(data[key], errors='coerce')
    valid = data.longitude.between(-180, 180) & data.latitude.between(-90, 90)
    data['valid_coordinates'] = valid
    data.to_csv(out / 'geography_coverage.csv', index=False)
    if not valid.any():
        return 'skipped: no verified numeric material coordinates; country centroids not substituted'
    fig, ax = plt.subplots(figsize=(8, 5))
    for i, (hid, group) in enumerate(data[valid].groupby('haplotype_id', sort=False)):
        ax.scatter(group.longitude, group.latitude, label=f'{hid} (n={len(group)})',
                   color=HAP_COLORS[i % len(HAP_COLORS)], s=45, alpha=.65)
    ax.set(xlabel='Longitude (degrees)', ylabel='Latitude (degrees)', title='Available material coordinates')
    ax.grid(alpha=.2)
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', frameon=False)
    save(fig, out / 'figures', 'geography')
    return f'ok: {int(valid.sum())} materials with coordinates; no coordinate imputation'


def canonical_value(value):
    value = str(value).strip()
    if value.upper() in {'', 'NA', 'N/A', 'NAN', 'NULL', '.', '?'}:
        return None
    try:
        from decimal import Decimal, InvalidOperation
        number = Decimal(value)
        if number.is_finite():
            return format(number.normalize(), 'f')
    except InvalidOperation:
        pass
    return value


def phenotype_summary(unique, data, ids):
    data = data.copy()
    data['value'] = data.value.map(canonical_value)
    data = data.dropna(subset=['value']).drop_duplicates(['sample_id', 'trait', 'value'])
    conflicts = data.groupby(['sample_id', 'trait']).value.nunique()
    bad = conflicts[conflicts > 1].reset_index()[['sample_id', 'trait']]
    marked = data.merge(bad.assign(conflict=True), on=['sample_id', 'trait'], how='left')
    clean = marked[marked.conflict.isna()].drop(columns='conflict')
    matched = clean.merge(unique[['sample_id', 'haplotype_id']], on='sample_id', how='inner', validate='many_to_one')
    rows = []
    for trait in sorted(data.trait.unique()):
        categories = sorted(data.loc[data.trait == trait, 'value'].unique())
        for hid in ids:
            subset = matched[(matched.trait == trait) & (matched.haplotype_id == hid)]
            n = len(subset)
            for category in categories:
                count = int((subset.value == category).sum())
                rows.append({'trait': trait, 'haplotype_id': hid, 'category': category,
                             'count': count, 'effective_n': n, 'percent': count / n * 100 if n else np.nan})
    return pd.DataFrame(rows, columns=['trait', 'haplotype_id', 'category', 'count', 'effective_n', 'percent']), bad, matched


def phenotypes(unique, path, ids, out):
    if path is None or unique is None:
        return 'skipped: phenotype file or unique material assignments unavailable'
    data = read(path, ['sample_id', 'trait', 'value'])
    summary, conflicts, matched = phenotype_summary(unique, data, ids)
    summary.to_csv(out / 'phenotype_percentages.csv', index=False)
    conflicts.to_csv(out / 'phenotype_conflicts_excluded.csv', index=False)
    matched.to_csv(out / 'phenotype_effective_records.csv', index=False)
    pd.DataFrame({'sample_id': sorted(set(data.sample_id) - set(unique.sample_id))}).to_csv(out / 'phenotype_unmatched_samples.csv', index=False)
    for number, (trait, group) in enumerate(summary.groupby('trait'), 1):
        fig, ax = plt.subplots(figsize=(max(7, len(ids) * .7 + 3), 4.8))
        bottom = np.zeros(len(ids))
        for i, category in enumerate(group.category.unique()):
            rows = group[group.category == category].set_index('haplotype_id').reindex(ids)
            values = rows.percent.to_numpy()
            ax.bar(range(len(ids)), values, bottom=bottom, label=category,
                   color=RESPONSE_COLORS.get(category, HAP_COLORS[i % len(HAP_COLORS)]))
            for j, row in enumerate(rows.itertuples()):
                if row.count and row.percent >= 8:
                    ax.text(j, bottom[j] + row.percent / 2, str(row.count), ha='center', va='center', fontsize=12)
            bottom += np.nan_to_num(values)
        ns = group.drop_duplicates('haplotype_id').set_index('haplotype_id').effective_n
        for i, hid in enumerate(ids):
            ax.text(i, 102, f'n={ns.get(hid, 0)}', ha='center', fontsize=12)
        ax.set_xticks(range(len(ids)), ids, rotation=45)
        ax.set(ylim=(0, 115), ylabel='Percent of materials with a valid record (%)', title=trait)
        ax.legend(title='Raw category', bbox_to_anchor=(1.02, 1), loc='upper left', frameon=False)
        save(fig, out / 'figures', f'phenotype_{number:02d}')
    return f'ok: descriptive categories; {len(matched)} valid material × trait records; {len(conflicts)} conflicts excluded'


def proteins(haps, sites, counts, assignments, manifest, input_dir, out, domains_path):
    if manifest.get('cds_phase_valid') is not True:
        return 'skipped: reference CDS phase/length not validated', None
    fasta = input_dir / 'reference_cds.fasta'
    if not fasta.exists():
        return 'skipped: reference CDS unavailable', None
    sequence = ''.join(x.strip() for x in fasta.read_text().splitlines() if not x.startswith('>'))
    try:
        sequences, reference = reconstruct_proteins(haps, sites, sequence, manifest['strand'])
    except ValueError as error:
        return f'skipped: {error}', None
    types, mapping = {}, {}
    for hid, peptide in sequences.items():
        if peptide not in types:
            types[peptide] = f'P{len(types) + 1}'
        mapping[hid] = types[peptide]
    pd.DataFrame([{'haplotype_id': h, 'protein_type': p} for h, p in mapping.items()]).to_csv(out / 'cds_to_protein.csv', index=False)
    (out / 'retained_site_proteins.fasta').write_text(''.join(f'>{pid}\n{seq}\n' for seq, pid in types.items()), encoding='utf-8')
    changed = [i for i in range(len(reference)) if any(seq[i] != reference[i] for seq in types)]
    aa = pd.DataFrame([{'haplotype_id': pid, **{str(i + 1): seq[i] for i in changed}} for seq, pid in types.items()])
    aa.to_csv(out / 'protein_variable_residues.csv', index=False)
    matrix_plot(aa, [str(i + 1) for i in changed], [str(i + 1) for i in changed], out / 'figures', 'protein_matrix', 'Protein types defined by retained SNPs', amino=True)
    collapsed = counts.assign(haplotype_id=counts.haplotype_id.map(mapping)).groupby(['haplotype_id', 'population'], as_index=False)['count'].sum()
    frequencies(collapsed, list(types.values()), out, 'protein')
    unique, excluded = materials(assignments.assign(haplotype_id=assignments.haplotype_id.map(mapping)) if assignments is not None else None)
    excluded.to_csv(out / 'protein_ambiguous_materials.csv', index=False)
    if unique is not None:
        unique.to_csv(out / 'protein_material_assignments.csv', index=False)
    fig, ax = plt.subplots(figsize=(9, 3.5))
    length = len(reference.rstrip('*'))
    ax.hlines(0, 1, length, color='#555555', linewidth=5)
    if domains_path:
        domains = read(domains_path, ['name', 'start', 'end', 'source'])
        for row in domains.itertuples():
            start, end = int(row.start), int(row.end)
            if not 1 <= start <= end <= length or not row.source.strip():
                raise ValueError('Domain needs a source and valid 1-based reference-protein coordinates')
            ax.add_patch(Rectangle((start, -.12), max(end - start, 1), .24, color='#7A9FC4'))
            ax.text((start + end) / 2, -.28, row.name, ha='center', fontsize=12)
        domains.to_csv(out / 'protein_domains.csv', index=False)
    for i in changed:
        ax.vlines(i + 1, 0, .5, color='#D55E00')
        ax.text(i + 1, .55, str(i + 1), rotation=60, fontsize=12)
    ax.set(xlim=(0, length + 1), ylim=(-.6, 1), xlabel='Reference amino acid position (1-based)',
           title='Variable residues from retained CDS SNPs', yticks=[])
    save(fig, out / 'figures', 'protein_positions')
    return f'ok: {len(types)} retained-site protein types; {len(changed)} variable residues', association(unique, 'protein')


def analyze(args):
    input_dir = Path(args.input)
    manifest = json.loads((input_dir / 'manifest.json').read_text(encoding='utf-8'))
    if manifest.get('role') not in ['target', 'tf'] or manifest.get('strand') not in ['+', '-'] or not manifest.get('assembly'):
        raise ValueError('Manifest needs role target|tf, strand +|-, and verified assembly')
    haps = read(input_dir / 'haplotypes.csv', ['haplotype_id'])
    sites = read(input_dir / 'sites.csv', ['site_id', 'position', 'alleles'])
    for key in sites.site_id:
        if key in haps:
            haps[key] = haps[key].str.strip().str.upper()
    counts = read(input_dir / 'counts.csv', ['haplotype_id', 'population', 'count'])
    assignments = read(input_dir / 'assignments.csv', ['sample_id', 'haplotype_id', 'population']) if (input_dir / 'assignments.csv').exists() else None
    validate(haps, sites, counts, assignments)
    keys = sites.site_id.tolist()
    if manifest.get('site_columns', keys) != keys:
        raise ValueError('Manifest site order differs from matrix definition')
    out = empty_output(args.out)
    setup()
    modules, tests = {}, []
    unique, excluded = materials(assignments)
    excluded.to_csv(out / 'ambiguous_materials.csv', index=False)
    if unique is not None:
        unique.to_csv(out / 'material_assignments.csv', index=False)
    if haps.empty or not keys:
        modules['core'] = 'skipped: no retained variants/types; no monomorphic cohort invented'
    else:
        haps.to_csv(out / 'haplotype_matrix.csv', index=False)
        sites.to_csv(out / 'retained_sites.csv', index=False)
        labels = sites.relative_position.tolist() if manifest['role'] == 'target' and 'relative_position' in sites else sites.position.tolist()
        matrix_plot(haps, keys, labels, out / 'figures', 'haplotype_matrix', 'Promoter haplotypes (relative to gene boundary)' if manifest['role'] == 'target' else 'CDS haplotypes (genomic positions)')
        table = frequencies(counts, haps.haplotype_id.tolist(), out, 'haplotype')
        modules['core'] = 'ok'
        modules['network'] = network(haps, keys, table, out)
        modules['geography'] = geography(unique, input_dir / 'metadata.csv', out)
        modules['phenotypes'] = phenotypes(unique, Path(args.phenotypes) if args.phenotypes else None, haps.haplotype_id.tolist(), out)
        tests.append(association(unique, 'promoter' if manifest['role'] == 'target' else 'CDS'))
        if manifest['role'] == 'tf':
            modules['protein'], test = proteins(haps, sites, counts, assignments, manifest, input_dir, out, Path(args.domains) if args.domains else None)
            if test:
                tests.append(test)
    adjust_bh(tests)
    pd.DataFrame(tests).to_csv(out / 'population_tests.csv', index=False)
    sources = list(input_dir.glob('*.csv')) + [input_dir / 'manifest.json'] + list(input_dir.glob('*.fasta'))
    sources += [Path(x) for x in [args.phenotypes, args.domains] if x]
    audit = {'created_utc': datetime.now(timezone.utc).isoformat(), 'manifest': manifest, 'modules': modules,
             'tests': tests, 'inputs': [{'file': p.name, 'sha256': sha256(p)} for p in sources]}
    write_json(out / 'run_audit.json', audit)
    report = [f'# {manifest.get("symbol", manifest.get("gene_id", "Gene"))} 单倍型分析', '',
              '**合成练习数据，不代表真实生物学结果。**' if manifest.get('synthetic') else '输入版本与基因身份见 run_audit.json。', '',
              f'- 分支：{"启动子" if manifest["role"] == "target" else "CDS／蛋白"}；参考版本：{manifest["assembly"]}',
              f'- 返回类型：{len(haps)}；保留位点：{len(keys)}；返回等位拷贝：{int(counts["count"].sum())}',
              f'- 唯一归属材料：{len(unique) if unique is not None else "未知"}；多类型材料：{len(excluded)}', '',
              '## 模块状态', '', *[f'- {k}: {v}' for k, v in modules.items()], '', '## 如何读图', '',
              '矩阵每行是一个类型，每列是一个保留位点。群体图分别回答“这个群体有哪些类型”和“这个类型来自哪些群体”；分母均为返回的等位拷贝。网络边数字是保留位点状态差异，不能直接解释为祖先、进化方向或选择。', '',
              '表型图只描述有效材料记录，n=0 是缺数据。类别名称按原始数据保留，不能自行赋予抗性强弱。群体检验只使用唯一归属的独立材料，P/q 不证明自然选择或基因功能，也未控制群体亲缘关系。', '',
              '## 参数与证据边界', '', json.dumps(manifest, ensure_ascii=False, indent=2), '',
              '图件及其源表位于当前目录。蛋白类型仅由保留 SNP 定义，未保留位点可能隐藏额外差异。缺少坐标、表型、序列或 motif 输入的模块应补齐数据后另行执行；未运行的扩展不能称为阴性结果。', '',
              '统计源表 population_tests.csv；BH 仅校正本次实际运行的群体检验，不是所有历史探索的联合校正。']
    (out / 'report_zh.md').write_text('\n'.join(report) + '\n', encoding='utf-8')
    print(f'Analysis complete: {out / "report_zh.md"}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True)
    parser.add_argument('--out', required=True)
    parser.add_argument('--phenotypes')
    parser.add_argument('--domains')
    try:
        analyze(parser.parse_args())
    except (ValueError, KeyError, FileNotFoundError) as error:
        print(f'Analysis stopped: {error}', file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
