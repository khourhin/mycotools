from pathlib import Path

import click
import pandas as pd
import requests


def attributes_to_json(attributes):
    return {x.split("=")[0]: x.split("=")[1] for x in attributes.split(";")}


def gff_to_df(gff):
    df = pd.read_csv(
        gff,
        sep="\t",
        comment="#",
        names=["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"],
    )
    # Extract attributes and transform them to columns
    attributes_df = df.attributes.apply(attributes_to_json)
    df = pd.concat([df, pd.json_normalize(attributes_df)], axis=1)

    return df


def download_annotations():
    print("Downloading functionnal annotation and H37rv/Erdman correspondance...")
    func_annot = requests.get(
        "https://mycobrowser.epfl.ch/releases/3/get_file?dir=txt&file=Mycobacterium_tuberculosis_H37Rv.txt"
    )
    corr_h37rv = requests.get(
        "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=gff3&id=379026087&conwithfeat=on&hide-sequence=on&hide-cdd=on"
    )

    open("mycobrowser_functional_annotation.tsv", "wb").write(func_annot.content)
    open("H37rv_Erdman_correspondance.gff", "wb").write(corr_h37rv.content)


def join_func_and_gff(func_txt, gff, table_to_annot, id_col):
    # Get functionnal annotation
    func_df = pd.read_csv(func_txt, sep=None, engine="python")

    # Read GFF for mapping Rv/Erdman
    gff = gff_to_df(gff)

    # Keep only information for mapping Erdman IDs to Rv IDs
    df = gff.query("type == 'CDS'").loc[:, ["locus_tag", "Note"]]
    df.Note = df.Note.str.replace("simliar to ", "", regex=False).str.replace(
        " of M. tuberculosis H37Rv", "", regex=False
    )
    df.rename(columns={"Note": "Locus"}, inplace=True)

    print(f"{df.isna().sum().Locus} Erdman genes does not have a correspondence to Rv IDs")

    # Read table to annotate (removing redundant info)
    df_to_annot = pd.read_csv(table_to_annot, sep=None, engine="python").drop(["Name"], axis=1)

    # Merge functional annotation with Erdman IDs using Erdman/Rv correspondence
    df = pd.merge(func_df, df, on="Locus")
    df = pd.merge(df_to_annot, df, left_on=id_col, right_on="locus_tag", how="left")

    df = df.loc[
        :,
        [
            "Id",
            "Name",
            "old_locus_tag",
            "Function",
            "Functional_Category",
            "Product",
            "baseMean",
            "log2FoldChange",
            "pvalue",
            "padj",
        ],
    ]

    df.drop_duplicates(inplace=True)

    outfile = Path(table_to_annot).with_suffix(".ann.xlsx")

    writer = pd.ExcelWriter(outfile, engine="xlsxwriter")
    df.to_excel(writer, sheet_name="Sheet1", index=False)
    writer.sheets["Sheet1"]

    workbook = writer.book
    worksheet = writer.sheets["Sheet1"]

    color_code = {
        "virulence, detoxification, adaptation": workbook.add_format({"bg_color": "#b2df8a"}),
        "information pathways": workbook.add_format({"bg_color": "#e31a1c"}),
        "cell wall and cell processes": workbook.add_format({"bg_color": "#33a02c"}),
        "stable RNAs": workbook.add_format({"bg_color": "#1f78b4"}),
        "insertion seqs and phages": workbook.add_format({"bg_color": "#a6cee3"}),
        "PE/PPE": workbook.add_format({"bg_color": "#6a3d9a"}),
        "intermediary metabolism and respiration": workbook.add_format({"bg_color": "#ff7f00"}),
        "unknown": workbook.add_format({"bg_color": "#808080"}),
        "regulatory proteins": workbook.add_format({"bg_color": "#cab2d6"}),
        "conserved hypotheticals": workbook.add_format({"bg_color": "#fb9a99"}),
        "lipid metabolism": workbook.add_format({"bg_color": "#b15928"}),
    }

    for k in color_code:
        worksheet.conditional_format(
            f"$A$1:$J${df.shape[0]}",
            {"type": "formula", "criteria": f'=INDIRECT("E"&ROW())="{k}"', "format": color_code[k]},
        )
    workbook.close()


@click.command()
@click.argument("id_col")
@click.argument("tables_to_annot", nargs=-1)
@click.option(
    "--func_txt",
    help="the file downloaded from mycobrowser with functionnal annotations for H37rv and using 'RvXXXX' IDs. Functional annotation can be downloaded as tab separated file from here: https://mycobrowser.epfl.ch/releases. Automactically downloaded by default.",
    default="mycobrowser_functional_annotation.tsv",
)
@click.option(
    "--gff",
    help="a gff file where correspondence can be found between 'ERDMAN_XXXX' IDs and H37Rv 'RvXXXX' IDs using the 'Similar to' field. Initially, this correspondence file was found here: https://www.ncbi.nlm.nih.gov/nuccore/AP012340. Automatically downloaded by default.",
    default="H37rv_Erdman_correspondance.gff",
)
def main(tables_to_annot, id_col, func_txt, gff):
    """
    Annotate DR results with pretty functional annotation from Mycobrowser for M. tuberculosis ERDMAN samples.

    ID_COL: the name of the column in the TABLES_TO_ANNOT file(s) where ERDMAN ids formated as 'ERDMAN_XXXX' ('ERDMAN_3272' for example) can be found.

    TABLE_TO_ANNOT: the file we want to add IDs and functionnal annotation to.
    Usually this table is comming from DR analysis from the platform. A column with 'ERDMAN_XXXX' IDs  should be present.
    """
    download_annotations()
    for t in tables_to_annot:
        join_func_and_gff(func_txt, gff, t, id_col)


if __name__ == "__main__":
    main()
