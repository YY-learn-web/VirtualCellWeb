# backend/app/algorithms/scfea_core.py

import os
import time
import warnings
import torch
import numpy as np
import pandas as pd
import magic
from torch.autograd import Variable
from tqdm import tqdm
import matplotlib.pyplot as plt

# 引入同目录下的依赖
from .ClassFlux import FLUX
from .DatasetFlux import MyDataset
from .util import pearsonr

# 定义默认超参数
DEFAULT_PARAMS = {
    'LEARN_RATE': 0.008,
    'LAMB_BA': 1.0,
    'LAMB_NG': 1.0,
    'LAMB_CELL': 1.0,
    'LAMB_MOD': 1e-2,
    'EPOCH': 100
}


def myLoss(m, c, lamb1, lamb2, lamb3, lamb4, geneScale=None, moduleScale=None):
    # balance constrain
    total1 = torch.sum(torch.pow(c, 2), dim=1)

    # non-negative constrain
    error = torch.abs(m) - m
    total2 = torch.sum(error, dim=1)

    # sample-wise variation constrain
    diff = torch.pow(torch.sum(m, dim=1) - geneScale, 2)
    if sum(diff > 0) == m.shape[0]:
        total3 = torch.pow(diff, 0.5)
    else:
        total3 = diff

    # module-wise variation constrain
    if lamb4 > 0:
        corr = torch.zeros(m.shape[0])
        for i in range(m.shape[0]):
            corr[i] = pearsonr(m[i, :], moduleScale[i, :])
        corr = torch.abs(corr)
        total4 = torch.ones(m.shape[0]) - corr
    else:
        total4 = torch.zeros(m.shape[0])

    loss = torch.sum(lamb1 * total1) + torch.sum(lamb2 * total2) + \
           torch.sum(lamb3 * total3) + torch.sum(lamb4 * total4)
    return loss, torch.sum(lamb1 * total1), torch.sum(lamb2 * total2), \
        torch.sum(lamb3 * total3), torch.sum(lamb4 * total4)


def run_scfea_training(
        input_file: str,
        output_dir: str,
        assets_dir: str,
        sc_imputation: bool = False,
        epochs: int = 100
):
    """
    执行 scFEA 训练的主函数
    :param input_file: 用户上传的 CSV 文件路径
    :param output_dir: 结果输出文件夹（通常是 output/job_id/）
    :param assets_dir: 存放 module_gene, cmMat 等参考文件的目录
    :param sc_imputation: 是否进行 MAGIC 插补
    :param epochs: 训练轮数
    """

    # 路径配置
    moduleGene_file = os.path.join(assets_dir, 'module_gene_m168.csv')
    cm_file = os.path.join(assets_dir, 'cmMat_c70_m168.csv')
    cName_file = os.path.join(assets_dir, 'cName_c70_m168.csv')

    # 确保输出目录存在
    os.makedirs(output_dir, exist_ok=True)

    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")

    # ================= Load Data =================
    print("Loading data...")
    geneExpr = pd.read_csv(input_file, index_col=0)
    #geneExpr = geneExpr.T
    geneExpr = geneExpr * 1.0

    if sc_imputation:
        magic_operator = magic.MAGIC()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            geneExpr = magic_operator.fit_transform(geneExpr)

    if geneExpr.max().max() > 50:
        geneExpr = (geneExpr + 1).apply(np.log2)

    geneExprSum = geneExpr.sum(axis=1)
    stand = geneExprSum.mean()
    geneExprScale = geneExprSum / stand
    geneExprScale = torch.FloatTensor(geneExprScale.values).to(device)

    BATCH_SIZE = geneExpr.shape[0]

    # Load Module Gene
    moduleGene = pd.read_csv(moduleGene_file, sep=',', index_col=0)
    moduleLen = np.array([moduleGene.iloc[i, :].notna().sum() for i in range(moduleGene.shape[0])])

    # Find Overlap Genes
    module_gene_all = []
    for i in range(moduleGene.shape[0]):
        for j in range(moduleGene.shape[1]):
            if not pd.isna(moduleGene.iloc[i, j]):
                module_gene_all.append(moduleGene.iloc[i, j])
    module_gene_all = set(module_gene_all)
    data_gene_all = set(geneExpr.columns)
    gene_overlap = list(data_gene_all.intersection(module_gene_all))
    gene_overlap.sort()

    # Load Stoichiometry Matrix
    cmMat = pd.read_csv(cm_file, sep=',', header=None).values
    cmMat = torch.FloatTensor(cmMat).to(device)

    cName = None
    if os.path.exists(cName_file):
        cName = pd.read_csv(cName_file, sep=',', header=0).columns

    # ================= Process Data =================
    print("Processing data...")
    geneExpr = geneExpr[gene_overlap]
    gene_names = geneExpr.columns
    cell_names = geneExpr.index.astype(str)
    n_modules = moduleGene.shape[0]
    n_genes = len(gene_names)
    n_cells = len(cell_names)
    n_comps = cmMat.shape[0]

    geneExprDf = pd.DataFrame(columns=['Module_Gene'] + list(cell_names))
    emptyNode = []

    for i in range(n_modules):
        genes = moduleGene.iloc[i, :].values.astype(str)
        genes = [g for g in genes if g != 'nan']
        if not genes:
            emptyNode.append(i)
            continue
        temp = geneExpr.copy()
        # 这里的逻辑是：只保留属于该 module 的基因，其他置0
        cols_to_zero = [g for g in gene_names if g not in genes]
        if cols_to_zero:
            temp.loc[:, cols_to_zero] = 0

        temp = temp.T
        temp['Module_Gene'] = ['%02d_%s' % (i, g) for g in gene_names]
        geneExprDf = pd.concat([geneExprDf, temp], ignore_index=True, sort=False)

    geneExprDf.index = geneExprDf['Module_Gene']
    geneExprDf.drop('Module_Gene', axis='columns', inplace=True)
    X = geneExprDf.values.T
    X = torch.FloatTensor(X).to(device)

    # Constraint of module variation
    df = geneExprDf.copy()
    df.index = [i.split('_')[0] for i in df.index]
    df.index = df.index.astype(int)
    module_scale = df.groupby(df.index).sum().T
    module_scale = torch.FloatTensor(module_scale.values / moduleLen)

    # ================= Train NN =================
    print("Training neural network...")
    torch.manual_seed(16)
    net = FLUX(X, n_modules, f_in=n_genes, f_out=1).to(device)
    optimizer = torch.optim.Adam(net.parameters(), lr=DEFAULT_PARAMS['LEARN_RATE'])

    dataSet = MyDataset(X, geneExprScale, module_scale)
    train_loader = torch.utils.data.DataLoader(dataset=dataSet, batch_size=BATCH_SIZE, shuffle=False)

    loss_history = {'total': [], 'balance': [], 'neg': [], 'cell': [], 'mod': []}

    net.train()
    for epoch in tqdm(range(epochs)):
        loss_val, l1, l2, l3, l4 = 0, 0, 0, 0, 0
        for i, (X_batch_data, X_scale_data, m_scale_data) in enumerate(train_loader):
            X_batch = Variable(X_batch_data.float().to(device))
            X_scale_batch = Variable(X_scale_data.float().to(device))
            m_scale_batch = Variable(m_scale_data.float().to(device))

            out_m, out_c = net(X_batch, n_modules, n_genes, n_comps, cmMat)
            loss_batch, loss1, loss2, loss3, loss4 = myLoss(
                out_m, out_c,
                DEFAULT_PARAMS['LAMB_BA'], DEFAULT_PARAMS['LAMB_NG'],
                DEFAULT_PARAMS['LAMB_CELL'], DEFAULT_PARAMS['LAMB_MOD'],
                geneScale=X_scale_batch, moduleScale=m_scale_batch
            )

            optimizer.zero_grad()
            loss_batch.backward()
            optimizer.step()

            loss_val += loss_batch.item()
            l1 += loss1.item()
            l2 += loss2.item()
            l3 += loss3.item()
            l4 += loss4.item()

        loss_history['total'].append(loss_val)

    # Save Loss Plot
    plt.figure()
    plt.plot(loss_history['total'], '--', label='total')
    plt.legend()
    plt.savefig(os.path.join(output_dir, 'loss_curve.png'))
    plt.close()

    # ================= Predict & Save =================
    print("Generating results...")
    flux_matrix = np.zeros((n_cells, n_modules), dtype='f')
    balance_matrix = np.zeros((n_cells, n_comps), dtype='f')

    test_loader = torch.utils.data.DataLoader(dataset=dataSet, batch_size=1, shuffle=False)
    net.eval()

    with torch.no_grad():
        for i, (X_batch_data, _, _) in enumerate(test_loader):
            X_batch = Variable(X_batch_data.float().to(device))
            out_m, out_c = net(X_batch, n_modules, n_genes, n_comps, cmMat)
            flux_matrix[i, :] = out_m.cpu().numpy()
            balance_matrix[i, :] = out_c.cpu().numpy()

    # Save Files
    flux_filename = os.path.join(output_dir, "predicted_flux.csv")
    balance_filename = os.path.join(output_dir, "predicted_balance.csv")

    setF = pd.DataFrame(flux_matrix, columns=moduleGene.index, index=geneExpr.index.tolist())
    setF.to_csv(flux_filename)

    setB = pd.DataFrame(balance_matrix, index=setF.index)
    if cName is not None:
        setB.columns = cName
    else:
        setB.columns = [f"Comp_{i}" for i in range(n_comps)]
    setB.to_csv(balance_filename)

    print(f"scFEA Training completed. Files saved to {output_dir}")
    return flux_filename, balance_filename