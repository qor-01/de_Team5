import pandas as pd
import itertools
from tqdm import tqdm

# 결과 저장
all_nodes = []
all_edges = []

# 염색체 리스트
chromosomes = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

for chrom in tqdm(chromosomes):
    try:
        # 1. 파일 불러오기
        df = pd.read_csv(f"KOVA_{chrom}_net_p.csv")

        # 2. 유전자 정보가 있는 행만
        df_gene = df[df['vep.worst_csq_by_gene.gene_symbol'].notnull()]

        # 3. 유전자별 요약 정보
        grouped = df_gene.groupby("vep.worst_csq_by_gene.gene_symbol").agg({
            "CADD": "mean",
            "clinvar_pathogenic": "max",
            "gwas_class": "max"
        }).reset_index()

        # 3.1 CADD 기준 상위 100개만 남기기
        grouped = grouped.sort_values(by="CADD", ascending=False).head(100)

        # 4. 컬럼 정리
        grouped.rename(columns={"vep.worst_csq_by_gene.gene_symbol": "gene_symbol"}, inplace=True)

        # 5. 노드 누적
        all_nodes.append(grouped)

        # 6. 엣지 생성
        for (i, row1), (j, row2) in itertools.combinations(grouped.iterrows(), 2):
            gene1 = row1['gene_symbol']
            gene2 = row2['gene_symbol']
            cadd_diff = abs(row1['CADD'] - row2['CADD'])

            # 조건: CADD 차이가 작을수록 유사하다고 판단
            if cadd_diff < 0.5:
                all_edges.append((gene1, gene2, cadd_diff))

    except Exception as e:
        print(f"{chrom} 처리 중 오류 발생: {e}")

# 7. 노드 합치고 중복 제거
nodes_df = pd.concat(all_nodes, ignore_index=True).drop_duplicates(subset="gene_symbol")

# 8. 엣지 정리 및 개수 제한
edges_df = pd.DataFrame(all_edges, columns=["source", "target", "cadd_diff"])
edges_df = edges_df.sort_values(by="cadd_diff").head(10000)  # 가장 유사한 10,000개

# 9. 저장
nodes_df.to_csv("nodes_kova_filtered.csv", index=False)
edges_df.to_csv("edges_kova_filtered.csv", index=False)
