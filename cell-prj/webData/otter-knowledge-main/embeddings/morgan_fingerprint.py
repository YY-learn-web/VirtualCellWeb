from logging import getLogger

import numpy as np
import torch
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem


class MorganFingerprint:
    def __init__(self, shape=2048, radius=2):
        self.shape = shape
        self.radius = radius
        self.logger = getLogger(__name__)

    @staticmethod
    def canonicalize(smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            return Chem.MolToSmiles(mol, isomericSmiles=True)
        else:
            return smiles

    def morgan_finger_print(self, smile) -> torch.Tensor:
        try:
            smile = self.canonicalize(smile)
            mol = Chem.MolFromSmiles(smile)
            if mol is None:
                raise ValueError(f"无效的SMILES格式: '{smile}'")
            features_vec = AllChem.GetMorganFingerprintAsBitVect(
                mol, self.radius, nBits=self.shape
            )
            features = np.zeros((1,))
            DataStructs.ConvertToNumpyArray(features_vec, features)
        except Exception as e:
            self.logger.warning(f"rdkit not found this smiles for morgan: {smile} convert to all 0 features")
            with open('morgan_error_smiles.txt', 'a') as f:
                f.write(smile + '\n')
            features = np.zeros((self.shape,))
            # 抛出错误
            raise ValueError(f"无效的SMILES格式: '{smile}'. 错误详情: {str(e)}")
        return torch.tensor(features, dtype=torch.float32)

    def get_embedding(self, smiles):
        return [self.morgan_finger_print(s) for s in smiles]
