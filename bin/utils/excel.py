from typing import Iterable, TypeAlias

import polars as pl
from xlsxwriter import Workbook
import os


Columns: TypeAlias = str | Iterable[str]
PathLike: TypeAlias = str | os.PathLike


def write_summarized_genomes_to_xlsx(
    df: pl.DataFrame,
    output_file: PathLike,
    group_by: Columns,
    sort_order_columns: Columns,
    extra_frames: Iterable|dict = tuple(),
):
    format_kw = {
        "autofit": True,
        }
    # turn all this into an xlsx
    with Workbook(output_file) as wb:
        if df is not None and not df.is_empty():
            for sheet, frame in df.group_by(group_by):
                frame = frame.sort(sort_order_columns)
                frame = frame.drop(group_by)
                frame.write_excel(
                    workbook=wb,
                    worksheet=sheet[0],
                    **format_kw
                )
        for extra_frame in extra_frames:
            if extra_frame is None:
                continue
            # If dictionary of names to frames, get the name
            if isinstance(extra_frames, dict):
                name = extra_frame
                extra_frame = extra_frames[extra_frame]
            # if tuple of frames, get the name from group_by column
            else:
                name = str(extra_frame[group_by][0])
            if extra_frame is not None and not extra_frame.is_empty():
                extra_frame.write_excel(
                    workbook=wb,
                    worksheet=name,
                    **format_kw
                )
