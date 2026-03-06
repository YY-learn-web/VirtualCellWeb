import os
import pandas as pd
import numpy as np
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
import umap
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

sns.set(style="whitegrid", font_scale=1.1)


# === 新增：用于独立生成自定义对比火山图的函数 ===
def generate_custom_diff_plot(output_dir: str, c1: int, c2: int):
    """读取已有的 flux 和 cluster 数据，快速生成新的火山图"""
    flux_file = os.path.join(output_dir, "predicted_flux.csv")
    cluster_file = os.path.join(output_dir, "cell_clusters.csv")

    if not os.path.exists(flux_file) or not os.path.exists(cluster_file):
        raise FileNotFoundError("找不到分析结果文件，请确保主任务已完成。")

    flux = pd.read_csv(flux_file, index_col=0)
    cluster_series = pd.read_csv(cluster_file, index_col=0).iloc[:, 0]

    # 统一列名格式
    flux.columns = [str(c) if str(c).startswith("M_") else ("M_" + str(c).replace("M", "")) for c in flux.columns]

    # 标准化 (复用之前的逻辑)
    flux_log = np.log1p(flux)
    flux_std = flux_log.std()
    with np.errstate(divide='ignore', invalid='ignore'):
        flux_z = (flux_log - flux_log.mean()) / flux_std
    flux_z = flux_z.fillna(0).replace([np.inf, -np.inf], 0)

    group1 = cluster_series[cluster_series == c1].index
    group2 = cluster_series[cluster_series == c2].index

    if len(group1) == 0 or len(group2) == 0:
        raise ValueError(f"所选的 Cluster {c1} 或 Cluster {c2} 没有细胞。")

    res = []
    for m in flux_z.columns:
        val1 = flux_z.loc[group1, m]
        val2 = flux_z.loc[group2, m]
        delta = val1.mean() - val2.mean()
        try:
            stat, p = mannwhitneyu(val1, val2, alternative='two-sided')
        except ValueError:
            p = 1.0
        res.append((m, delta, p))

    diff_df = pd.DataFrame(res, columns=["module", "delta_flux", "pvalue"])
    diff_df["pvalue"] = diff_df["pvalue"].fillna(1.0)
    diff_df["FDR"] = multipletests(diff_df["pvalue"], method="fdr_bh")[1]

    # 画图
    diff_df["-log10FDR"] = -np.log10(diff_df["FDR"] + 1e-10)
    plt.figure(figsize=(6, 5))
    sns.scatterplot(data=diff_df, x="delta_flux", y="-log10FDR", s=20)
    plt.axhline(-np.log10(0.05), ls="--", c="red")
    plt.axvline(0, ls="--", c="grey")
    plt.xlabel(f"ΔFlux (Cluster {c1} - Cluster {c2})")
    plt.ylabel("-log10(FDR)")
    plt.title(f"Differential Flux Volcano (C{c1} vs C{c2})")
    plt.tight_layout()

    img_name = f"Differential_Flux_Volcano_{c1}_vs_{c2}.png"
    plt.savefig(os.path.join(output_dir, img_name), dpi=300)
    plt.close()

    csv_name = f"diff_flux_cluster{c1}_vs_{c2}.csv"
    diff_df.sort_values("FDR").to_csv(os.path.join(output_dir, csv_name), index=False)

    return img_name, csv_name


# === 修改：主分析函数接收 n_clusters 参数 ===
def run_downstream_analysis(flux_file: str, balance_file: str, output_dir: str, assets_dir: str, n_clusters: int = 4):
    print(f"Starting analysis with {n_clusters} clusters...")
    N_PCA = 10
    RANDOM_STATE = 0

    anno_file = os.path.join(assets_dir, "Human_M168_information.symbols.csv")
    anno = pd.read_csv(anno_file) if os.path.exists(anno_file) else pd.DataFrame()

    module_to_pathway = {
        'Glycolysis': ['M_1', 'M_2', 'M_3', 'M_4', 'M_5', 'M_6'],
        'TCA_Cycle': ['M_7', 'M_8', 'M_9', 'M_10', 'M_11', 'M_12', 'M_13', 'M_14'],
        'Pentose_Phosphate_Pathway': ['M_33'],
        'Glycogen_Metabolism': ['M_111'],
        'Fatty_Acid_Biosynthesis': ['M_34'],
        'Fatty_Acid_Oxidation': ['M_35'],
        'Steroid_Cholesterol_Biosynthesis': ['M_167', 'M_168', 'M_169'],
        'Serine_Glycine_Metabolism': [f'M_{i}' for i in range(15, 33)],
        'Aspartate_Metabolism': ['M_36', 'M_37', 'M_38', 'M_39', 'M_40'],
        'Alanine_Aspartate_Metabolism': ['M_41', 'M_42', 'M_43', 'M_44', 'M_45'],
        'Glutamine_Glutamate_Metabolism': ['M_48', 'M_49', 'M_50', 'M_51', 'M_52'],
        'BCAA_Metabolism': [f'M_{i}' for i in range(53, 61)],
        'Arginine_Proline_Metabolism': [f'M_{i}' for i in range(61, 69)],
        'Polyamine_Metabolism': ['M_69', 'M_70'],
        'Propanoate_Metabolism': ['M_46', 'M_47'],
        'Purine_Metabolism': [f'M_{i}' for i in range(133, 150)] + ['M_170'],
        'Pyrimidine_Metabolism': [f'M_{i}' for i in range(150, 167)] + ['M_171'],
        'Amino_Sugar_Metabolism': ['M_106', 'M_107', 'M_108', 'M_109', 'M_110'],
        'N_Glycan_Biosynthesis': [f'M_{i}' for i in range(112, 125)],
        'O_Glycan_Biosynthesis': ['M_125', 'M_126', 'M_127', 'M_128'],
        'Chondroitin_Dermatan_Biosynthesis': ['M_129', 'M_130', 'M_131'],
        'Heparan_Sulfate_Biosynthesis': ['M_132'],
        'Transporters': ['M_71', 'M_72', 'M_73', 'M_74', 'M_75', 'M_76', 'M_77', 'M_78', 'M_79', 'M_80', 'M_81', 'M_83',
                         'M_84', 'M_85', 'M_87', 'M_88', 'M_89', 'M_90', 'M_91', 'M_92', 'M_93', 'M_94', 'M_95', 'M_96',
                         'M_97', 'M_98', 'M_99', 'M_100', 'M_101', 'M_102', 'M_103', 'M_105']
    }
    module_to_pathway_flat = {m: p for p, ms in module_to_pathway.items() for m in ms}

    flux = pd.read_csv(flux_file, index_col=0)
    balance = pd.read_csv(balance_file, index_col=0)
    common_cells = flux.index.intersection(balance.index)
    flux = flux.loc[common_cells]
    balance = balance.loc[common_cells]

    flux.columns = [str(c) if str(c).startswith("M_") else ("M_" + str(c).replace("M", "")) for c in flux.columns]

    flux_log = np.log1p(flux)
    flux_std = flux_log.std()
    with np.errstate(divide='ignore', invalid='ignore'):
        flux_z = (flux_log - flux_log.mean()) / flux_std
    flux_z = flux_z.fillna(0).replace([np.inf, -np.inf], 0)

    pca = PCA(n_components=min(N_PCA, flux_z.shape[0], flux_z.shape[1]), random_state=RANDOM_STATE)
    X_pca = pca.fit_transform(flux_z)
    pca_df = pd.DataFrame(X_pca, index=flux_z.index, columns=[f"PC{i + 1}" for i in range(X_pca.shape[1])])
    pca_df.to_csv(os.path.join(output_dir, "PCA_coordinates.csv"))

    reducer = umap.UMAP(random_state=RANDOM_STATE)
    X_umap = reducer.fit_transform(flux_z)
    umap_df = pd.DataFrame(X_umap, index=flux_z.index, columns=["UMAP1", "UMAP2"])
    umap_df.to_csv(os.path.join(output_dir, "UMAP_coordinates.csv"))

    # === 使用传入的 n_clusters ===
    kmeans = KMeans(n_clusters=n_clusters, random_state=RANDOM_STATE, n_init=20)
    clusters = kmeans.fit_predict(flux_z)
    cluster_series = pd.Series(clusters, index=flux_z.index, name="cluster")
    cluster_series.to_csv(os.path.join(output_dir, "cell_clusters.csv"))

    flux_z["cluster"] = cluster_series
    flux_log["cluster"] = cluster_series
    cluster_flux_mean = flux_log.groupby("cluster").mean()

    # 默认生成 C0 vs C1 的火山图作为初始图
    default_volcano_img = "Differential_Flux_Volcano_0_vs_1.png"
    if n_clusters >= 2:
        default_volcano_img, _ = generate_custom_diff_plot(output_dir, 0, 1)

    module_pathway_series = pd.Series(module_to_pathway_flat)
    module_pathway_series = module_pathway_series.loc[
        module_pathway_series.index.intersection(flux_z.drop(columns="cluster").columns)]

    pathway_flux = flux_z.drop(columns="cluster").groupby(module_pathway_series, axis=1).mean()
    pathway_cluster_mean = pathway_flux.join(cluster_series).groupby("cluster").mean()

    stress_score = balance.abs().sum(axis=1)
    stress_score.name = "stress_score"
    flux_stress_corr = flux_log.drop(columns="cluster").corrwith(stress_score).sort_values(ascending=False)

    cell_meta = pd.concat([pca_df, umap_df, cluster_series, stress_score], axis=1)

    # 照搬 plot_figure.py 的标签映射逻辑
    module_display_name = {}
    if not anno.empty:
        for _, row in anno.iterrows():
            module = str(row.iloc[0])  # 强制取第一列，不检查列名

            # 标准化 M_x
            if module.startswith("M") and not module.startswith("M_"):
                module = "M_" + module[1:]
            elif module.isdigit():
                module = "M_" + module

            in_cpd = row.get("Compound_IN_name", "?")
            out_cpd = row.get("Compound_OUT_name", "?")
            pathway = module_to_pathway_flat.get(module, "Unknown")

            # 拼接成: Pathway | IN -> OUT
            module_display_name[module] = f"{pathway} | {in_cpd} -> {out_cpd}"

    # 1. UMAP by Cluster
    plt.figure(figsize=(6, 5))
    sns.scatterplot(data=cell_meta, x="UMAP1", y="UMAP2", hue="cluster", palette="tab10", s=12)
    plt.savefig(os.path.join(output_dir, "UMAP_by_cluster.png"), dpi=300)
    plt.close()

    # 2. UMAP by Stress
    plt.figure(figsize=(6, 5))
    plt.scatter(cell_meta["UMAP1"], cell_meta["UMAP2"], c=cell_meta["stress_score"], cmap="viridis", s=12)
    plt.colorbar(label="Metabolic stress score")
    plt.savefig(os.path.join(output_dir, "UMAP_by_stress.png"), dpi=300)
    plt.close()

    # 3. Cluster x Module Heatmap (完全还原你的调整)
    cluster_flux_named = cluster_flux_mean.rename(columns=module_display_name)
    plt.figure(figsize=(max(12, 0.35 * cluster_flux_named.shape[1]), 8))

    ax = sns.heatmap(
        cluster_flux_named,
        cmap="RdBu_r",
        center=0,
        xticklabels=True,
        yticklabels=True
    )
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontsize=6)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=9)
    plt.title("Cluster-wise metabolic module activity", fontsize=12)

    # 关键：照搬你原来的边距调整，防止长标签被截断
    plt.subplots_adjust(
        left=0.08,
        right=0.98,
        bottom=0.35,
        top=0.92
    )
    plt.savefig(os.path.join(output_dir, "Cluster_Module_Heatmap.png"), dpi=300)
    plt.close()

    # 4. Cluster x Pathway Heatmap
    plt.figure(figsize=(10, 6))
    sns.heatmap(pathway_cluster_mean, cmap="RdBu_r", center=0)
    plt.xlabel("Pathway")
    plt.ylabel("Cluster")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "Cluster_Pathway_Heatmap.png"), dpi=300)
    plt.close()

    # 5. Flux-Stress Correlation Barplot (还原标签)
    flux_stress_corr_named = flux_stress_corr.rename(index=module_display_name)
    top_corr = flux_stress_corr_named.head(20)
    plt.figure(figsize=(10, max(6, 0.4 * top_corr.shape[0])))
    ax = top_corr.plot(kind="bar", legend=False)
    ax.set_ylabel("Correlation with metabolic stress")
    plt.xticks(rotation=90, fontsize=7)
    plt.subplots_adjust(bottom=0.45)  # 同样增加底部边距防截断
    plt.savefig(os.path.join(output_dir, "Flux_Stress_Top20_vertical.png"), dpi=300)
    plt.close()

    return {
        "files": {
            "flux": flux_file,
            "pca": os.path.join(output_dir, "PCA_coordinates.csv"),
            "cluster": os.path.join(output_dir, "cell_clusters.csv")
        },
        "images": {
            "umap_cluster": "UMAP_by_cluster.png",
            "umap_stress": "UMAP_by_stress.png",
            "heatmap_module": "Cluster_Module_Heatmap.png",
            "heatmap_pathway": "Cluster_Pathway_Heatmap.png",
            "volcano": default_volcano_img,
            "stress_corr": "Flux_Stress_Top20_vertical.png"
        }
    }