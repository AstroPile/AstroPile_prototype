import sys
sys.path.append('../')
import torch
from pytorch_lightning import LightningDataModule
from torch.utils.data import Dataset, DataLoader
from utils import split_dataset, compute_dataset_statistics, normalize_sample
from typing import Any

class PhotoZWrapper(LightningDataModule):
    def __init__(self, 
                 dataset: Dataset, 
                 batch_size: int = 128, 
                 num_workers: int = 8, 
                 test_size: float = 0.2, 
                 split_method: str = 'naive', 
                 loading: str = 'iterated', 
                 feature_flag: str = 'image',
                 label_flag: str = 'z',
                 dynamic_range: bool = True):
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

        self.dataset = dataset
        self.batch_size = batch_size
        self.num_workers = num_workers
        self.test_size = test_size
        self.split_method = split_method
        self.loading = loading
        self.dynamic_range = dynamic_range
        self.feature_flag = feature_flag
        self.label_flag = label_flag

    def prepare_data(self):
        # Split the dataset
        self.train_dataset, self.test_dataset = split_dataset(self.dataset, split=self.split_method)

        # Compute the dataset statistics
        self.img_mean, self.img_std = compute_dataset_statistics(self.train_dataset, flag=self.feature_flag, loading=self.loading)
        self.z_mean, self.z_std = compute_dataset_statistics(self.train_dataset, flag=self.label_flag, loading=self.loading)

        train_test_split = self.train_dataset.train_test_split(test_size=self.test_size)
        self.train_dataset = train_test_split['train']
        self.val_dataset = train_test_split['test']

    def setup(self, stage=None):
        pass

    def collate_fn(self, batch):
        batch = torch.utils.data.default_collate(batch)
        batch['image'] = normalize_sample(batch['image'], self.img_mean, self.img_std, dynamic_range=self.dynamic_range).permute(0, 3, 1, 2) # dynamic range compression and z-score normalization
        batch['redshift'] = normalize_sample(batch['redshift'], self.z_mean, self.z_std, dynamic_range=False)
        return batch['image'], batch['redshift']

    def train_dataloader(self):
        return DataLoader(self.train_dataset, batch_size=self.batch_size, num_workers=self.num_workers, collate_fn=self.collate_fn)
    
    def val_dataloader(self):
        return DataLoader(self.val_dataset, batch_size=self.batch_size, num_workers=self.num_workers, collate_fn=self.collate_fn)
    
    def test_dataloader(self):
        return DataLoader(self.test_dataset, batch_size=self.batch_size, num_workers=self.num_workers, collate_fn=self.collate_fn)