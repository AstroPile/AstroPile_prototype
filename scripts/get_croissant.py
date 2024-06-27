# Retrieve Croissant description from Hugging Face Hub
import json
import os
import requests
from typing import List

API_TOKEN = os.environ["HG_API_TOKEN"]  # Get API Token


def query(url: str):
    headers = {"Authorization": f"Bearer {API_TOKEN}"}
    response = requests.get(url, headers=headers)
    return response.json()


def get_datasets() -> List[str]:
    dataset_urls = "https://huggingface.co/api/datasets/?author=AstroPile"
    data = query(dataset_urls)
    available_datasets = []
    for dataset in data:
        if "croissant" in dataset["tags"]:
            available_datasets.append(dataset["id"])
    return available_datasets


def get_croissant(dataset_id: str):
    dataset_name = dataset_id.split("/")[1]
    croissant_url = f"https://huggingface.co/api/datasets/{dataset_id}/croissant"
    croissant_data = query(croissant_url)
    with open(f"{dataset_name}.json", "w") as output_file:
        json.dump(croissant_data, output_file)


if __name__ == "__main__":
    for dataset_id in get_datasets():
        print(f"Get Croissant for {dataset_id}")
        get_croissant(dataset_id)
