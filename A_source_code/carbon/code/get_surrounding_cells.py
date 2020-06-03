# ******************************************************
## Copyright 2019, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

def do(icell, dum_asc, type="icell"):
    '''
    Returns list of neighboring cells
    '''
    nbs = list()
    irow,icol = dum_asc.get_row_col_from_index(icell)
    for row in range(irow-1, irow+2):
        for col in range(icol-1, icol+2):
            if not (col==icol and row==irow):
                if type=="icell":
                  nbs.append(dum_asc.get_index_from_row_col(row, col))
                elif type=="irow_icol":
                  nbs.append([row, col])
    return nbs