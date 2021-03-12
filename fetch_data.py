"""
Scripts to fech data from the web for my projects. Depends on having activated the GSheets API.
Info on the website of gspread: https://gspread.readthedocs.io/en/latest/oauth2.html
"""


def gsheet_to_csv(gsheet_name, worksheet, csv_filename):
    import os

    from gspread import service_account
    from pandas import DataFrame

    """
    Get as a csv data from a google sheet.
    
    Parameters
    ----------
    csv_filename: path and filename of the csv.
    gsheet_name: name of the spreadsheet.
    worksheet: name of the worksheet.
      
    Returns
    -------
    """
    if "/" in csv_filename:
        assert os.path.exists(
            os.path.dirname(csv_filename)
        ), f"Directory {os.path.dirname(csv_filename)} does not exist."

    gc = service_account()
    sh = gc.open(gsheet_name)
    worksheet = sh.worksheet(worksheet)
    DataFrame(worksheet.get_all_records()).to_csv(csv_filename)
    return


def gsheet_to_csv_url(gsheet_name, worksheet, csv_filename):
    import os

    from gspread import service_account
    from pandas import DataFrame

    """
    Get as a csv data from a google sheet.
    
    Parameters
    ----------
    csv_filename: path and filename of the csv.
    gsheet_name: name of the spreadsheet.
    worksheet: name of the worksheet.
      
    Returns
    -------
    """
    if "/" in csv_filename:
        assert os.path.exists(
            os.path.dirname(csv_filename)
        ), f"Directory {os.path.dirname(csv_filename)} does not exist."

    gc = service_account()
    sh = gc.open_by_url(gsheet_name)
    worksheet = sh.worksheet(worksheet)
    DataFrame(worksheet.get_all_records()).to_csv(csv_filename)
    return


def df_to_data_sheet(df, gsheet_name, worksheet):
    """
    Send  dataframe df to worksheet in gsheet_name google sheet.

    Parameters
    ----------
    def: dataframe  to send to gsheets.
    gsheet_name: name of the spreadsheet.
    worksheet: name of the worksheet.

    Returns
    -------
    """
    from gspread import service_account
    from pandas import DataFrame

    assert isinstance(df, DataFrame), "df is not a pandas DataFrame"
    assert isinstance(gsheet_name, str), "gsheet_name is not a string"
    assert isinstance(worksheet, str), "gsheet_name is not a string"
    gc = service_account()
    sh = gc.open(gsheet_name)
    worksheet = sh.worksheet(worksheet)
    worksheet.update([df.columns.values.tolist()] + df.values.tolist())
    return
