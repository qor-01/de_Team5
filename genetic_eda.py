import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import numpy as np
import os

# 한글 폰트 설정
font_path = 'C:/Windows/Fonts/malgun.ttf'
fontprop = fm.FontProperties(fname=font_path)
plt.rcParams['font.family'] = fm.FontProperties(fname=font_path).get_name()

# 공통 시각화 함수
def plot_top_bar(df, col, palette='viridis', top_n=10):
    top_vals = df[col].value_counts().nlargest(top_n)
    plt.figure(figsize=(10, 4))
    sns.barplot(x=top_vals.values, y=top_vals.index, palette=palette)
    plt.title(f'Top {top_n} frequent values in {col}')
    plt.xlabel('Count')
    plt.ylabel(col)
    plt.tight_layout()
    plt.show()

def plot_pie_chart(series, title):
    plt.figure(figsize=(6, 6))
    plt.pie(series, labels=series.index, autopct='%1.1f%%', startangle=140)
    plt.title(title)
    plt.tight_layout()
    plt.show()

def plot_histogram(series, title, color):
    plt.figure(figsize=(8, 4))
    sns.histplot(series, bins=30, kde=True, color=color)
    plt.title(title)
    plt.xlabel("Length")
    plt.ylabel("Count")
    plt.tight_layout()
    plt.show()

# 분석할 컬럼
columns = [
'kova_AC', 'kova_AF', 'jpn_AC', 'jpn_AF', 'gnomad_eas_AF',
    'gnomad_eur_AF', 'gnomad_all_AF', 'chimpanzee_AF', 'clinvar_CLNSIG',
    'clinvar_pathogenic', 'gwas_pvalue', 'gwas_class', 'gtex_qval',
    'gtex_class', 'CADD', 'FunSeq2', 'ReMM', 'LINSIGHT', 'sift',
    'phyloP100way', 'phastCons100way',
    "vep.most_severe_consequence", "vep.worst_consequence_term",
    "vep.worst_csq_by_gene.impact", "vep.worst_csq_by_gene.lof",
    "vep.worst_csq_by_gene.gene_symbol", "vep.worst_csq_by_gene.gene_id",
    "vep.worst_csq_by_gene.hgnc_id", "vep.worst_csq_by_gene.amino_acids",
    "vep.worst_csq_by_gene.hgvsc", "vep.worst_csq_by_gene.hgvsp"
]

# 크로모좀 목록
chromosomes = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

# 루프 실행
for chrom in chromosomes:
    file_path = f"1_KOVA.v7.{chrom}.csv"
    if not os.path.exists(file_path):
        print(f"{file_path} not found. Skipping.")
        continue

    print(f"\n🔍 Processing {file_path} ...")
    df = pd.read_csv(file_path, low_memory=False)

    df = df[[col for col in columns if col in df.columns]]

    for col in ["vep.most_severe_consequence", "vep.worst_consequence_term",
                "vep.worst_csq_by_gene.impact", "vep.worst_csq_by_gene.lof"]:
        if col in df.columns:
            plot_top_bar(df, col)

    if "vep.worst_csq_by_gene.gene_symbol" in df.columns:
        top_genes = df["vep.worst_csq_by_gene.gene_symbol"].value_counts().nlargest(10)
        plot_pie_chart(top_genes, "Top 10 Most Frequent Gene Symbols")

    if "vep.worst_csq_by_gene.amino_acids" in df.columns:
        df["aa_length"] = df["vep.worst_csq_by_gene.amino_acids"].astype(str).apply(len)
        plot_histogram(df["aa_length"], "Distribution of Amino Acid Change Length", "orange")

    for col, color in [("vep.worst_csq_by_gene.hgvsc", "blue"), ("vep.worst_csq_by_gene.hgvsp", "green")]:
        if col in df.columns:
            df[f"{col}_length"] = df[col].astype(str).apply(len)
            plot_histogram(df[f"{col}_length"], f"Distribution of {col.split('.')[-1].upper()} Length", color)

    if "vep.most_severe_consequence" in df.columns and "vep.worst_csq_by_gene.impact" in df.columns:
        plt.figure(figsize=(12, 6))
        sns.countplot(data=df, x="vep.most_severe_consequence", hue="vep.worst_csq_by_gene.impact", palette="Set2")
        plt.title("Most Severe Consequence by VEP Impact")
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.show()

    # 결측치, 고유값, 기타 정보 출력
    print(df.head(3))
    print(df.isnull().sum())

    # veplist 고유값 확인
    veplist = columns
    for i in veplist:
        if i in df.columns:
            num_i = df[i].unique()
            print(f"{i}: {len(num_i)} unique values")
            print(num_i)

    # 집단별 AF 시각화
    af_columns = ["kova_AF", "jpn_AF", "gnomad_eas_AF", "gnomad_eur_AF", "gnomad_all_AF", "chimpanzee_AF"]
    af_present = [col for col in af_columns if col in df.columns]

    if af_present:
        df_melted = df[af_present].melt(var_name="Population", value_name="Allele Frequency")
        plt.figure(figsize=(10, 6))
        sns.boxplot(x="Population", y="Allele Frequency", data=df_melted)
        plt.xlabel("집단", fontproperties=fontprop)
        plt.ylabel("대립유전자 빈도 (AF)", fontproperties=fontprop)
        plt.title("집단별 대립유전자 빈도 분포 비교", fontproperties=fontprop)
        plt.xticks(rotation=30)
        plt.show()

    # 병원성 여부에 따른 대립유전자 빈도
    if "clinvar_pathogenic" in df.columns and "kova_AF" in df.columns:
        plt.figure(figsize=(8, 6))
        sns.boxplot(x=df["clinvar_pathogenic"], y=df["kova_AF"])
        plt.xlabel("병원성 여부 (0: 비병원성, 1: 병원성)", fontproperties=fontprop)
        plt.ylabel("대립유전자 빈도 (kova_AF)", fontproperties=fontprop)
        plt.title("병원성 여부에 따른 대립유전자 빈도 비교", fontproperties=fontprop)
        plt.show()

    # GWAS 시각화
    if "gwas_pvalue" in df.columns:
        df = df[df["gwas_pvalue"] > 0]
        df["log_gwas_pvalue"] = -np.log10(df["gwas_pvalue"])

        plt.figure(figsize=(8, 5))
        sns.scatterplot(x=df["kova_AF"], y=df["log_gwas_pvalue"], alpha=0.5)
        plt.xlabel("대립유전자 빈도 (kova_AF)", fontproperties=fontprop)
        plt.ylabel("-log10(GWAS P-value)", fontproperties=fontprop)
        plt.title("GWAS P-value와 대립유전자 빈도 관계", fontproperties=fontprop)
        plt.show()

    if "phyloP100way" in df.columns and "phastCons100way" in df.columns:
        plt.figure(figsize=(8, 5))
        sns.scatterplot(x=df["phyloP100way"], y=df["phastCons100way"], alpha=0.5)
        plt.xlabel("phyloP100way 보존 점수", fontproperties=fontprop)
        plt.ylabel("phastCons100way 보존 점수", fontproperties=fontprop)
        plt.title("보존 점수 비교 (phyloP vs phastCons)", fontproperties=fontprop)
        plt.show()

    if "gwas_class" in df.columns and "log_gwas_pvalue" in df.columns:
        plt.figure(figsize=(10, 6))
        sns.boxplot(x="gwas_class", y="log_gwas_pvalue", data=df)
        plt.xlabel("GWAS 질병 연관성 클래스", fontproperties=fontprop)
        plt.ylabel("-log10(GWAS P-value)", fontproperties=fontprop)
        plt.title("GWAS 클래스별 P-value 분포 비교", fontproperties=fontprop)
        plt.xticks(rotation=30, fontproperties=fontprop)
        plt.show()

        plt.figure(figsize=(10, 6))
        sns.violinplot(x="gwas_class", y="log_gwas_pvalue", data=df)
        plt.xlabel("GWAS 질병 연관성 클래스", fontproperties=fontprop)
        plt.ylabel("-log10(GWAS P-value)", fontproperties=fontprop)
        plt.title("GWAS 클래스별 P-value 분포 (Violin Plot)", fontproperties=fontprop)
        plt.xticks(rotation=30, fontproperties=fontprop)
        plt.show()
