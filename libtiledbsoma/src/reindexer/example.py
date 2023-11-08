import cellxgene_census
import tiledbsoma as soma
import os
from time import perf_counter
#soma.pytiledbsoma.config_logging("debug")
print(os.getpid())

def main():
    census_S3_latest = dict(census_version="latest")
    census_local_copy = dict(uri="/Users/brobatmili/projects/census_data/")
    t1 = perf_counter()
    with cellxgene_census.open_soma(**census_local_copy) as census:
        with census["census_data"]["homo_sapiens"].axis_query(
                measurement_name="RNA", obs_query=soma.AxisQuery(
                    value_filter="""tissue_general == 'eye'""")

        ) as query:
            query.to_anndata(X_name="raw")
    t2 = perf_counter()
    print(f"End to end time {t2 - t1}")

main()
