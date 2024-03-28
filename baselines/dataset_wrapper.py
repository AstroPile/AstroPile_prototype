import sys
sys.path.append('../')
import torch
from pytorch_lightning import LightningDataModule
from torch.utils.data import Dataset, DataLoader
from dataset_utils import split_dataset, compute_dataset_statistics, normalize_sample, get_nested
from typing import Any

class DatasetWrapper(LightningDataModule):
    def __init__(self, 
                 train_dataset: Dataset, 
                 test_dataset: Dataset,
                 feature_flag: str = 'image',
                 label_flag: str = 'z',
                 feature_dynamic_range: bool = True,
                 feature_z_score: bool = True,
                 label_dynamic_range: bool = False,
                 label_z_score: bool = True,
                 batch_size: int = 128, 
                 num_workers: int = 8, 
                 test_size: float = 0.2, 
                 loading: str = 'iterated', 
                 ):
        super().__init__()
        """
        Initializes the data module with a dataset that is already loaded, setting up parameters
        for data processing and batch loading.

        Parameters:
        - dataset (Dataset): The pre-loaded dataset, expected to be a torch.utils.data.Dataset with images in size B x C x H x W.
        - batch_size (int): The size of each data batch for loading.
        - num_workers (int): Number of subprocesses to use for data loading.
        - test_size (float): The proportion of the dataset to reserve for testing.
        - split_method (str): Strategy for splitting the dataset ('naive' implemented).
        - loading (str): Approach for loading the dataset ('full' or 'iterated').
        - feature_flag (str): The key in the dataset corresponding to the image data.
        - label_flag (str): The key in the dataset corresponding to the redshift data.
        - dynamic_range (bool): Flag indicating whether dynamic range compression should be applied.
        """

        self.train_dataset = train_dataset
        self.test_dataset = test_dataset

        self.feature_flag = feature_flag
        self.feature_dynamic_range = feature_dynamic_range
        self.feature_z_score = feature_z_score

        self.label_flag = label_flag
        self.label_dynamic_range = label_dynamic_range
        self.label_z_score = label_z_score

        self.batch_size = batch_size
        self.num_workers = num_workers
        self.test_size = test_size
        self.loading = loading

        self.prepare_data_per_node = False

    def prepare_data(self):
        # Compute the dataset statistics
        if self.feature_z_score:
            self.feature_mean, self.feature_std = compute_dataset_statistics(self.train_dataset, flag=self.feature_flag, loading=self.loading)
        if self.label_z_score:
            self.label_mean, self.label_std = compute_dataset_statistics(self.train_dataset, flag=self.label_flag, loading=self.loading)

        # Split the dataset
        train_test_split = self.train_dataset.train_test_split(test_size=self.test_size)
        self.train_dataset = train_test_split['train']
        self.val_dataset = train_test_split['test']

    def setup(self, stage=None):
        pass

    def collate_fn(self, batch):
        batch = torch.utils.data.default_collate(batch)
        x = normalize_sample(get_nested(batch, self.feature_flag), self.feature_mean, self.feature_std, dynamic_range=self.feature_dynamic_range, z_score=self.feature_z_score) # dynamic range compression and z-score normalization
        y = normalize_sample(get_nested(batch, self.label_flag), self.feature_mean, self.feature_std, dynamic_range=self.label_dynamic_range, z_score=self.label_z_score)
        return x, y

    def train_dataloader(self):
        return DataLoader(self.train_dataset, batch_size=self.batch_size, num_workers=self.num_workers, collate_fn=self.collate_fn)
    
    def val_dataloader(self):
        return DataLoader(self.val_dataset, batch_size=self.batch_size, num_workers=self.num_workers, collate_fn=self.collate_fn)
    
    def test_dataloader(self):
        return DataLoader(self.test_dataset, batch_size=self.batch_size, num_workers=self.num_workers, collate_fn=self.collate_fn)
    
    
