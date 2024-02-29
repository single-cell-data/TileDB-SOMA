from time import perf_counter

import cellxgene_census

import tiledbsoma as soma

census_S3_latest = dict(census_version="2023-10-23")


def main():
    t1 = perf_counter()
    with cellxgene_census.open_soma(**census_S3_latest) as census:
        with census["census_data"]["homo_sapiens"].axis_query(
            measurement_name="RNA",
            obs_query=soma.AxisQuery(value_filter="""tissue_general == 'eye'"""),
        ) as query:
            query.to_anndata(X_name="raw")
    t2 = perf_counter()
    print(f"End to end time {t2 - t1}")


main()
