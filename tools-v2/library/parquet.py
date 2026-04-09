import os
import shutil
import polars as pl


class BufferedParquetWriter:

    def __init__(self, dir_name, n_rows=3, mode="w"):

        self.dir_name = dir_name
        self.n_rows = n_rows
        self.dataframe = None

        match mode:
            case "w":
                if os.path.exists(self.dir_name):
                    shutil.rmtree(self.dir_name)
                os.makedirs(self.dir_name, exist_ok=True)
                self.partition_num = 0
                self.schema = None
            case "a":
                if os.path.exists(self.dir_name):
                    self.partition_num = max(int(i.split(".")[0]) for i in os.listdir(dir_name))
                    self.schema = pl.scan_parquet(self.dir_name).collect_schema()
                else:
                    os.makedirs(self.dir_name, exist_ok=True)
                    self.partition_num = 0
                    self.schema = None
            case _:
                raise Exception()

    def get_partition_name(self):

        self.partition_num += 1
        return f"{self.dir_name}/{self.partition_num:0>8}.parquet"
    
    def add_dataframe(self, dataframe):

        if self.schema is None:
            self.schema = dataframe.schema
        if self.dataframe is None:
            self.dataframe = dataframe
        else:
            assert dataframe.schema == self.schema
            self.dataframe = pl.concat([self.dataframe, dataframe])
        
        *partitions, self.dataframe = self.dataframe.iter_slices(self.n_rows)
        for partition in partitions:
            partition.write_parquet(self.get_partition_name())

    def dump_buffer(self):

        if self.dataframe.shape[0] != 0:
            self.dataframe.write_parquet(self.get_partition_name())
            self.dataframe = self.schema.to_frame()

def repartition_parquets(input_dir, output_dir, n_rows=10_000):

    pl.scan_parquet(input_dir).collect_schema()
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir, exist_ok=True)
    
    partition_num = 0
    partitions = sorted(f"{input_dir}/{i}" for i in os.listdir(input_dir))
    for n, partition in enumerate(partitions):
        if n == 0:
            dataframe = pl.read_parquet(partition)
        else:
            dataframe = pl.concat([dataframe, pl.read_parquet(partition)])
        *partitions, dataframe = dataframe.iter_slices(n_rows)
        for partition in partitions:
            partition.write_parquet(f"{output_dir}/{partition_num:0>8}.parquet")
            partition_num += 1
    if dataframe.shape[0] != 0:
        dataframe.write_parquet(f"{output_dir}/{partition_num:0>8}.parquet")
