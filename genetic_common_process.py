import pandas as pd
import numpy as np
df = pd.read_csv(f"KOVA_chrY.csv",low_memory=False)
pd.set_option("display.max_columns", None)

selected_columns = [
    "kova_AF", "kova_AC", "jpn_AF", "jpn_AC", "gnomad_eas_AF", "gnomad_eas_AC",
    "gnomad_eur_AF", "gnomad_eur_AC", "gnomad_all_AF", "gnomad_all_AC", "chimpanzee_AF", "chimpanzee_AC",
    "CADD", "FunSeq2", "ReMM", "LINSIGHT", "sift",
    "clinvar_CLNSIG", "clinvar_pathogenic",
    "gwas_pvalue", "gwas_class",
    "phyloP100way", "phastCons100way",
    "gtex_qval", "gtex_class"
    , "vep.most_severe_consequence", "vep.worst_consequence_term",
    "vep.worst_csq_by_gene.impact", "vep.worst_csq_by_gene.lof",
    "vep.worst_csq_by_gene.gene_symbol","vep.worst_csq_by_gene.gene_id",
    "vep.worst_csq_by_gene.hgnc_id","vep.worst_csq_by_gene.amino_acids",
    "vep.worst_csq_by_gene.hgvsc","vep.worst_csq_by_gene.hgvsp"
]

df = df[selected_columns]

# 1. AF/AC 계열 (집단별 변이 빈도 및 개수) → 결측값 0으로 대체
af_ac_cols = [
        'jpn_AF', 'jpn_AC',
        'gnomad_eas_AF', 'gnomad_eas_AC',
        'gnomad_eur_AF', 'gnomad_eur_AC',
        'gnomad_all_AF', 'gnomad_all_AC',
        'chimpanzee_AF', 'chimpanzee_AC'
    ]
df[af_ac_cols] = df[af_ac_cols].fillna(0)

    # 2. 기능 예측 점수 → 결측값 유지 (분석 요소로 사용)
functional_scores = ['CADD', 'FunSeq2', 'ReMM', 'LINSIGHT']
    # 유지: 따로 처리 안 함

    # 3. sift → 데이터 편향으로 인해 Unknown 또는 NaN 유지
df['sift'] = df['sift'].fillna('Unknown')

    # 4. clinvar 계열 → Unknown으로 대체
df['clinvar_CLNSIG'] = df['clinvar_CLNSIG'].fillna('Unknown')
df['clinvar_pathogenic'] = df['clinvar_pathogenic'].fillna('Unknown')

    # 5. gwas_pvalue → 결측값 유지 또는 Unknown (여기서는 NaN 유지)
df['gwas_pvalue'] = df['gwas_pvalue'].fillna('Unknown')  # 선택사항

    # 6. gwas_class → 최빈값으로 대체
if df['gwas_class'].notna().any():
    mode_value = df['gwas_class'].mode()[0]
    df['gwas_class'] = df['gwas_class'].fillna(mode_value)

    # 7. 진화적 보존 점수 → 결측값 0으로 대체
df['phyloP100way'] = df['phyloP100way'].fillna(0)
df['phastCons100way'] = df['phastCons100way'].fillna(0)

    # 8. GTEx 관련 → gtex_qval은 0으로 대체, gtex_class는 삭제
df['gtex_qval'] = df['gtex_qval'].fillna(0)
df['gtex_class'] = df['gtex_class'].fillna('Unknwon')

    # 9. VEP 관련
    # - vep.worst_consequence_term: Unknown으로 대체
df['vep.worst_consequence_term'] = df['vep.worst_consequence_term'].fillna('Unknown')

    # - vep.worst_csq_by_gene.impact, lof: Unknown으로 대체
df['vep.worst_csq_by_gene.impact'] = df['vep.worst_csq_by_gene.impact'].fillna('Unknown')
df['vep.worst_csq_by_gene.lof'] = df['vep.worst_csq_by_gene.lof'].fillna('Unknown')

    # - vep.worst_csq_by_gene.* : 결측값 → Unknown
vep_cols = [
        'vep.worst_csq_by_gene.gene_symbol',
        'vep.worst_csq_by_gene.gene_id',
        'vep.worst_csq_by_gene.hgnc_id',
        'vep.worst_csq_by_gene.amino_acids',
        'vep.worst_csq_by_gene.hgvsc',
        'vep.worst_csq_by_gene.hgvsp'
    ]
df[vep_cols] = df[vep_cols].fillna('Unknown')

    # (필요 시, vep.most_severe_consequence 는 이미 결측 없음)
df.to_csv('KOVA_chrY_p.csv', index=False, encoding='utf-8')