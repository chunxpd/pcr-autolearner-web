import os
import yaml


class Config:
    def __init__(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        with open(f"{cur_dir}/../config.yaml", "r", encoding="utf-8") as handle:
            self.config = yaml.load(handle, Loader=yaml.FullLoader)
