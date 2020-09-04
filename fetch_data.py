def gsheet_to_csv(gsheet_name,worksheet,csv_filename):
    from gspread import service_account
    from pandas import DataFrame
    import os 
    '''
    Get as a csv data from a google sheet.
    
    Parameters
    ----------
    csv_filename: path and filename of the csv.
    gsheet_name: name of the spreadsheet.
    worksheet: name of the worksheet.
      
    Returns
    -------
    '''
    if '/' in csv_filename:
        assert os.path.exists(os.path.dirname(csv_filename)),f'Directory {os.path.dirname(csv_filename)} does not exist.'
    
    gc = service_account()
    sh = gc.open(gsheet_name)
    worksheet=sh.worksheet(worksheet)
    DataFrame(worksheet.get_all_records()).to_csv(csv_filename)
    return
