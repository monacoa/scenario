from pyxll import xlcAlert
import sys
import datetime
from win32com.client import constants as const
from Tkinter import *
from sc_elab.core.SwpCurve import dict_segm2, Segm, Curve, BootstrappedCurve




def popup_messagebox(msg):
    xlcAlert(msg)
#----
def drawBox(xla, spessore , rTopLeft = 0, cTopLeft = 0,rBottomRight=0,cBottomRight=0, Colore=0):
    if (rTopLeft <= 0) or(cTopLeft <= 0) or (rBottomRight <= 0)or (cBottomRight <= 0):
        msg = "Le coordinate del box devono essere maggiori di zero"
        popup_messagebox(msg)
        sys.exit()

    RR = xla.Range(xla.Cells(rTopLeft, cTopLeft), xla.Cells(rBottomRight, cBottomRight))
    RR.Borders(const.xlDiagonalDown).LineStyle = const.xlNone
    RR.Borders(const.xlDiagonalUp).LineStyle = const.xlNone
    RR.Borders(const.xlEdgeLeft).LineStyle = const.xlContinuous
    RR.Borders(const.xlEdgeLeft).Weight = spessore
    RR.Borders(const.xlEdgeLeft).ColorIndex = Colore

    RR.Borders(const.xlEdgeTop).LineStyle = const.xlContinuous
    RR.Borders(const.xlEdgeTop).Weight = spessore
    RR.Borders(const.xlEdgeTop).ColorIndex = Colore

    RR.Borders(const.xlEdgeBottom).LineStyle = const.xlContinuous
    RR.Borders(const.xlEdgeBottom).Weight = spessore
    RR.Borders(const.xlEdgeBottom).ColorIndex = Colore

    RR.Borders(const.xlEdgeRight).LineStyle = const.xlContinuous
    RR.Borders(const.xlEdgeRight).Weight = spessore
    RR.Borders(const.xlEdgeRight).ColorIndex = Colore


def drawLine (xla, rTopLeft, cTopLeft, rBottomRight, cBottomRight, hor, spessore):
    (xla.Range(xla.Cells(rTopLeft, cTopLeft), xla.Cells(rBottomRight, cBottomRight))).Select()
    if (hor =="o"):
        xla.Selection.Borders(const.xlEdgeBottom).LineStyle = const.xlContinuous
        xla.Selection.Borders(const.xlEdgeBottom).Weight = spessore
        xla.Selection.Borders(const.xlEdgeBottom).ColorIndex = const.xlAutomatic

    else:
        xla.Selection.Borders(const.xlEdgeRight).LineStyle = const.xlContinuous
        xla.Selection.Borders(const.xlEdgeRight).Weight = spessore
        xla.Selection.Borders(const.xlEdgeRight).ColorIndex = const.xlAutomatic
#----------
def formatTestataCurva(xla,nRiga,nColonna,nLarghezza,testo):
    (xla.Range(xla.Cells(nRiga, nColonna), xla.Cells(nRiga, nColonna + nLarghezza - 1))).Select()
    xla.Selection.HorizontalAlignment = const.xlCenter
    xla.Selection.VerticalAlignment = const.xlBottom
    xla.Selection.WrapText = False
    xla.Selection.Orientation = 0
    xla.Selection.AddIndent = False
    xla.Selection.IndentLevel = 0
    xla.Selection.ShrinkToFit = False
    xla.Selection.ReadingOrder = const.xlContext
    xla.Selection.MergeCells = False
    xla.Selection.Merge()
    xla.Selection.Font.ColorIndex = 2
    xla.Selection.Font.Bold = True
    xla.Selection.Interior.ColorIndex = 55
    xla.Selection.Interior.Pattern = const.xlSolid
    xla.Selection.Value = testo
#----------

def readCurvesNames(xla, s, rangeStart, direzione, distanza):
    r = xla.Range(rangeStart)
    curveL = []

    if direzione.lower() == "o":
        i = 0
        while (r.Value != None):
            i += 1

            nomeCurva = r.Value
            curveL.append((nomeCurva,i))

            nCols = r.Columns.Count
            row   = r.Row
            col   = r.Column
            r = xla.Range(xla.Cells(row, col + distanza), xla.Cells(row, col + (nCols - 1) + distanza))
    else:
        j = 1
        nomeCurva = r.Value
        curveL.append((nomeCurva, j))

        while (r.Value != None):
            nCols = r.Columns.Count
            row   = r.Row
            col   = r.Column
            r = xla.Range(xla.Cells(row + 1, col), xla.Cells(row + 1, col))
            # porto avanti ancora per controllare

            if (r.Value == None):
                r = xla.Range(xla.Cells(row  + distanza, col), xla.Cells(row + distanza, col))
                if r.Value != None:
                    j += 1
                    nomeCurva = r.Value
                    curveL.append((nomeCurva, j))

    return curveL

def findRigthPlaceBootCurveSeg(xla, r, distCurve, dir="O"):
    rOut = None
    if dir == "v" :
        if (r.Value == None): return r
        nCols = r.Columns.Count
        row   = r.Row
        col   = r.Column
        j = 0
        while (r.Value != None):
            j += 1
            r = xla.Range(xla.Cells(row + j, col), xla.Cells(row + j, col))
            #porto avanti ancora per controllare
            if (r.Value == None):  r = xla.Range(xla.Cells(row + j + distCurve , col), xla.Cells(row + j+distCurve, col))
        r = xla.Range(xla.Cells(row + j + distCurve, col), xla.Cells(row+j+distCurve, col))
    else:
        while (r.Value != None):

            nCols = r.Columns.Count
            row =r.Row
            col =r.Column
            r = xla.Range(xla.Cells(row, col + distCurve), xla.Cells(row, col + (nCols-1) + distCurve))
    rOut = r
    if (rOut == None):
        msg = "Unable to compute the output range for your curve"
        print msg
        sys.exit()

    return rOut