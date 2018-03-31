#!/usr/bin/env python

"""

"""

__author__ = "Ruben Acuna"
__copyright_ = "Copyright(c) 2011, ASU iGEM Team"

################################### IMPORTS ####################################
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio import Entrez
from Bio.Seq import Seq
import pickle
import os
import sys
import wx

from Bio import Restriction

import ArrayGeneration
import Config
import Updater
from utilities import UtilBioBrick
from utilities import UtilGeneral

################################################################################

class TextRedirector():
    def __init__(self, textctrlReceiver):
        self.textctrlReceiver = textctrlReceiver

    def write(self, text):
        self.textctrlReceiver.WriteText(text)
        
class DialogSpacerGenSettings_PAM(wx.Dialog):
    def __init__(self, parent):
        wx.Dialog.__init__(self, parent, title="Spacer Generation Settings: PAMs", size=(355, 400))
        self.settings_pickle = list(pickle.load(open(Config.SETTINGS)))

        #labels for lists
        self.statictextRESites = wx.StaticText(self, label="Available PAMs", pos=(14, 5))

        #build list of PAMs
        self.listctrlSelectedPAMs = wx.ListCtrl(self, style=wx.LC_REPORT, size=(225, 140), pos=(14, 25))
        self.listctrlSelectedPAMs.InsertColumn(0, 'Name')
        self.listctrlSelectedPAMs.InsertColumn(1, 'Sequence')
        self.listctrlSelectedPAMs.InsertColumn(2, 'Position')
        self.listctrlSelectedPAMs.SetColumnWidth(1, 90)
        self.listctrlSelectedPAMs.SetColumnWidth(2, 45)
        
        for i in self.settings_pickle[1][1][1:][::-1]:
            list_init = self.listctrlSelectedPAMs.InsertStringItem(0, str(i[1][0]))
            self.listctrlSelectedPAMs.SetStringItem(list_init, 1, str(i[1][1]))
            if i[1][2] == 1:
                self.listctrlSelectedPAMs.SetStringItem(list_init, 2, "5'")
            else:
                self.listctrlSelectedPAMs.SetStringItem(list_init, 2, "3'")
                
        self.y_coord = 184
        
        #text field for name
        self.statictextPAMName = wx.StaticText(self, label="Name:", pos=(14, self.y_coord))
        self.textctrlNameInput = wx.TextCtrl(self, size=(100, 20), pos=(84, self.y_coord - 2))
        self.textctrlNameInput.SetMaxLength(20)

        self.y_coord += 25
        
        #text field for entering sequence
        self.statictextPAMSeq = wx.StaticText(self, label="Sequence:", pos=(14, self.y_coord))
        self.textctrlPAMInput = wx.TextCtrl(self, size=(100, 20), pos=(84, self.y_coord - 2))

        self.y_coord += 25
        
        #radio buttons for position
        self.statictextPAMSeq = wx.StaticText(self, label="Position:", pos=(14, self.y_coord))

        self.radio_PAM_l = wx.RadioButton(self, -1, "5'", pos=(84, self.y_coord + 2), style=wx.RB_GROUP)
        self.radio_PAM_r = wx.RadioButton(self, -1, "3'", pos=(114, self.y_coord + 2))

        self.y_coord += 25
        
        #add / remove buttons
        self.buttonAddPAM = wx.Button(self, label="Add PAM", pos=(14, self.y_coord), size=(159, 25))
        self.buttonRemovePAM = wx.Button(self, label="Remove PAM", pos=(179, self.y_coord), size=(159, 25))

        self.Bind(wx.EVT_BUTTON, self.add_PAM, self.buttonAddPAM)
        self.Bind(wx.EVT_BUTTON, self.remove_PAM, self.buttonRemovePAM)

        self.y_coord += 30
        
        #dropdown
        self.statictextdropdown = wx.StaticText(self, label="PAM to use:", pos=(14, self.y_coord))
        dropdown_l = [a[1][0] for a in self.settings_pickle[1][1]]
        
        self.comboPAM = wx.ComboBox(self, -1, pos=(84, self.y_coord), size=(150, -1), choices=dropdown_l, style=wx.CB_READONLY)
        self.comboPAM.SetSelection(self.settings_pickle[1][0])

        self.y_coord += 25
        
        #save button
        self.buttonSave = wx.Button(self, label="Save settings", pos=(83, self.y_coord), size=(159, 25))
        self.Bind(wx.EVT_BUTTON, self.save_settings, self.buttonSave)

    def update_PAMs(self):

        self.listctrlSelectedPAMs.ClearAll()
        self.listctrlSelectedPAMs.InsertColumn(0, 'Name')
        self.listctrlSelectedPAMs.InsertColumn(1, 'Sequence')
        self.listctrlSelectedPAMs.InsertColumn(2, 'Position')
        self.listctrlSelectedPAMs.SetColumnWidth(1, 90)
        self.listctrlSelectedPAMs.SetColumnWidth(2, 45)

        for i in self.settings_pickle[1][1][1:][::-1]:
            list_init = self.listctrlSelectedPAMs.InsertStringItem(0, str(i[1][0]))
            self.listctrlSelectedPAMs.SetStringItem(list_init, 1, str(i[1][1]))
            if i[1][2] == 1:
                self.listctrlSelectedPAMs.SetStringItem(list_init, 2, "5'")
            else:
                self.listctrlSelectedPAMs.SetStringItem(list_init, 2, "3'")            
            
        self.comboPAM.Clear()
        dropdown_l = [a[1][0] for a in self.settings_pickle[1][1]]
        self.comboPAM.AppendItems(dropdown_l)
        self.comboPAM.SetSelection(self.settings_pickle[1][0])

    def add_PAM(self, event):
        valid_nucleotides = set(['A', 'C', 'T', 'G', 'N'])
        
        pam_name = self.textctrlNameInput.GetValue()
        pam_seq = self.textctrlPAMInput.GetValue().upper()
        pam_pos = self.radio_PAM_l.GetValue() 

        if pam_name != '':
            if pam_seq != '':
                if set(pam_seq) - valid_nucleotides == set([]):
                    PAM = [1, [pam_name, pam_seq, pam_pos]]
                else:
                    print 'Invalid nucleotides in input'
                    return
            else:
                print 'Enter a non-blank name and sequence'
                return
        else:
            print 'Enter a non-blank name and sequence'
            return
        self.settings_pickle[1][1].append(PAM)
        self.update_PAMs()

    def remove_PAM(self, event):
        site = self.listctrlSelectedPAMs.GetFocusedItem() + 1

        if site != -1:
            self.settings_pickle[1][1] = self.settings_pickle[1][1][:site] + self.settings_pickle[1][1][site + 1:]
            self.update_PAMs()
            
    def save_settings(self, event):
        self.settings_pickle[1][0] = self.comboPAM.GetCurrentSelection()
        wfile = open(Config.SETTINGS, "w")
        pickle.dump(self.settings_pickle, wfile)
        wfile.close()
        print "Settings saved"
        
class DialogSpacerGenSettings_RE(wx.Dialog):
    def __init__(self, parent):
        wx.Dialog.__init__(self, parent, title="Spacer Generation Settings: Restriction", size=(355, 315))
        self.settings_pickle = list(pickle.load(open(Config.SETTINGS)))

        self.radio_RE_yes = wx.RadioButton(self, -1, 'Filter restriction sites', pos=(10, 5), style=wx.RB_GROUP)
        self.radio_RE_no = wx.RadioButton(self, -1, "Don't filter", pos=(150, 5))
        self.radio_RE_no.SetValue(self.settings_pickle[0][0])
        
        #labels for lists
        self.statictextRESites = wx.StaticText(self, label="Selected", pos=(14, 25))
        self.statictextRESites = wx.StaticText(self, label="Available", pos=(180, 25))

        #build list of currently selected sites
        self.listctrlSelectedSites = wx.ListCtrl(self, style=wx.LC_REPORT, size=(157, 140), pos=(14, 45))
        self.update_sites()

        #build list of available sites
        self.listctrlAvailableSites = wx.ListCtrl(self, style=wx.LC_REPORT, size=(157, 140), pos=(180, 45))
        self.listctrlAvailableSites.InsertColumn(0, 'Name')
        self.listctrlAvailableSites.InsertColumn(1, 'Sequence')
        self.listctrlAvailableSites.SetColumnWidth(0, 45)
        self.listctrlAvailableSites.SetColumnWidth(1, 90)

        self.availableREs = [getattr(Restriction, x) for x in dir(Restriction) if type(getattr(Restriction, x)).__name__ == "RestrictionType"]

        for i in xrange(len(self.availableREs)):
            self.listctrlAvailableSites.InsertStringItem(i, str(self.availableREs[i]))
            self.listctrlAvailableSites.SetStringItem(i, 1, self.availableREs[i].site)

        #text field for entering re by name
        self.statictextRESites = wx.StaticText(self, label="Pick by Name: ", pos=(14, 194))
        self.textctrlSiteInput = wx.TextCtrl(self, size=(253, 20), pos=(84, 192))
        self.textctrlSiteInput.SetMaxLength(10)

        #buttons for adding and removing res
        self.buttonAddSite = wx.Button(self, label="Add Site", pos=(13,217), size=(159, 25))
        self.buttonRemoveSite = wx.Button(self, label="Remove Site", pos=(179, 217), size=(159, 25))

        self.Bind(wx.EVT_BUTTON, self.add_site, self.buttonAddSite)
        self.Bind(wx.EVT_BUTTON, self.remove_site, self.buttonRemoveSite)

        #save button
        self.buttonSave = wx.Button(self, label="Save settings", pos=(83, 252), size=(159, 25))
        self.Bind(wx.EVT_BUTTON, self.save_settings, self.buttonSave)
                
    def update_sites(self):

        self.listctrlSelectedSites.ClearAll()
        self.listctrlSelectedSites.InsertColumn(0, 'Name')
        self.listctrlSelectedSites.InsertColumn(1, 'Sequence')
        self.listctrlSelectedSites.SetColumnWidth(0, 45)
        self.listctrlSelectedSites.SetColumnWidth(1, 90)
        
        for RE in self.settings_pickle[0][1][::-1]:
            pos = self.listctrlSelectedSites.InsertStringItem(0, str(RE))
            self.listctrlSelectedSites.SetStringItem(pos, 1, RE.site)

    def add_site(self, event):
        if self.textctrlSiteInput.GetValue() != "":
            site = self.textctrlSiteInput.GetValue()
        elif self.listctrlAvailableSites.GetFocusedItem() != -1:
            site = str(self.availableREs[self.listctrlAvailableSites.GetFocusedItem()])
        else:
            site = ""
            
        try:
            self.settings_pickle[0][1].append(getattr(Restriction, site))
            self.update_sites()
            
        except AttributeError:
            print "Invalid site"

    def remove_site(self, event):
        site = self.listctrlSelectedSites.GetFocusedItem()

        if site != -1:
            self.settings_pickle[0][1] = self.settings_pickle[0][1][:site] + self.settings_pickle[0][1][site + 1:]
            self.update_sites()
            
    def save_settings(self, event):
        self.settings_pickle[0][0] = self.radio_RE_no.GetValue()

        wfile = open(Config.SETTINGS, "w")
        pickle.dump(self.settings_pickle, wfile)
        wfile.close()
        print "Settings saved"

class DialogSpacersFinalize(wx.Dialog):
    def __init__(self, parent, spacers, repeat):
        wx.Dialog.__init__(self, parent, title="Finalize Spacers", size=(320, 300))

        self.spacers = spacers
        self.repeat = repeat

        self.statictextConsole = wx.StaticText(self, label="Available Spacers:", pos=(10, 5))
        self.listSpacers = wx.ListBox(parent=self, id=wx.ID_ANY, size=(294, 100), pos=(10, 20))
        self.listSpacers.Items = [str(spacer) for spacer in spacers]
        self.Bind(wx.EVT_LISTBOX, self.OnClickListSpacers, self.listSpacers)

        self.statictextLayout = wx.StaticText(self, label="Sequence Layout:", pos=(10, 130))
        self.radioFormRS = wx.RadioButton(self, -1, 'RS', pos=(10, 150), style=wx.RB_GROUP)
        self.radioFormRSR = wx.RadioButton(self, -1, 'RSR', pos=(10, 170))
        self.Bind(wx.EVT_RADIOBUTTON, self.OnSelectRadioLayout, self.radioFormRS)
        self.Bind(wx.EVT_RADIOBUTTON, self.OnSelectRadioLayout, self.radioFormRSR)

        self.statictextLayout = wx.StaticText(self, label="Sequence Package:", pos=(210, 130))
        self.radioPackageNone = wx.RadioButton(self, -1, 'None', pos=(210, 150), style=wx.RB_GROUP)
        self.radioPackageBB = wx.RadioButton(self, -1, 'BioBrick', pos=(210, 170))
        self.Bind(wx.EVT_RADIOBUTTON, self.OnSelectRadioPackage, self.radioPackageNone)
        self.Bind(wx.EVT_RADIOBUTTON, self.OnSelectRadioPackage, self.radioPackageBB)

        self.statictextSpacer = wx.StaticText(self, label="Result:", pos=(10, 220))
        self.textctrlSpacer = wx.TextCtrl(self, size=(295, 20), pos=(10, 240))
        self.textctrlSpacer.SetFont(wx.Font(pointSize=10, family=wx.FONTFAMILY_MODERN, style=wx.NORMAL, weight=wx.FONTWEIGHT_NORMAL))
        self.textctrlSpacer.SetEditable(False)
        self.textctrlSpacer.SetForegroundColour((156,156,156)) #srsly?

    def UpdateResult(self):

        if self.listSpacers.GetSelection() == -1:
            self.textctrlSpacer.SetValue("")
            return

        sequence = self.listSpacers.Items[self.listSpacers.GetSelection()]

        if self.radioFormRS.GetValue():
            sequence = self.repeat + sequence
        else: #RSR
            sequence = self.repeat + sequence + self.repeat

        if self.radioPackageBB.GetValue():
            sequence = UtilBioBrick.makeBioBrick(sequence)

        self.textctrlSpacer.SetValue(str(sequence.upper()))

    def OnClickListSpacers(self, event):
        self.UpdateResult()

    def OnSelectRadioLayout(self, event):
        self.UpdateResult()

    def OnSelectRadioPackage(self, event):
        self.UpdateResult()

class PanelBrowser(wx.Panel):

    def BindGnomes(self):
        #bind outside GUI
        #parent.Bind(wx.EVT_LISTBOX, self.OnClickListGenomes, listGenomes)
        self.parent.listGenomes.Bind(wx.EVT_LISTBOX, self.OnClickListGenomes)

    def __init__(self, parent, dimensions, location, listGenomes, crisprIndex, dictCRISPRs):
        wx.Panel.__init__(self, parent, size=dimensions, pos=location)

        self.parent = parent
        self.BindGnomes()
        self.listGenome = listGenomes
        self.crisprIndex = crisprIndex
        self.dictCRISPRs = dictCRISPRs

        self.casGenes = pickle.load(open("data/casgeneDB.txt"))

        self.staticboxRepeat = wx.StaticBox(self, -1, "Organism", size=(379, 180-2), pos=(5, 10))

        self.statictextName = wx.StaticText(self, label="Name:", pos=(14, 30+3))
        self.textctrlName = wx.TextCtrl(self, size=(317, 20), pos=(60, 30))
        self.textctrlName.SetEditable(False)
        self.textctrlName.SetForegroundColour((156, 156, 156)) #srsly?

        self.statictextRefSeq = wx.StaticText(self, label="RefSeq:", pos=(14, 55+3))
        self.textctrlRefSeq = wx.TextCtrl(self, size=(317, 20), pos=(60, 55))
        self.textctrlRefSeq.SetEditable(False)
        self.textctrlRefSeq.SetForegroundColour((156, 156, 156)) #srsly?

        self.statictextTaxID = wx.StaticText(self, label="TaxID:", pos=(14, 80+3))
        self.textctrlTaxID = wx.TextCtrl(self, size=(317, 20), pos=(60, 80))
        self.textctrlTaxID.SetEditable(False)
        self.textctrlTaxID.SetForegroundColour((156, 156, 156)) #srsly?

        self.statictextAccessions = wx.StaticText(self, label="GenBank Accessions:", pos=(14, 105+3))
        self.listAccessions = wx.ListBox(self, id=wx.ID_ANY, size=(363, 55), pos=(14, 125))

        #array information
        self.staticboxRepeat = wx.StaticBox(self, -1, "Arrays", size=(379, 180-17), pos=(5, 195))
        self.statictextArAccessions = wx.StaticText(self, label="Accessions:", pos=(14, 215))
        self.listArrays = wx.ListBox(self, id=wx.ID_ANY, size=(100, 118), pos=(14, 230+2))
        parent.Bind(wx.EVT_LISTBOX, self.OnClickArray, self.listArrays)

        self.statictextRepeat = wx.StaticText(self, label="Repeat:", pos=(125, 215))
        self.textctrlRepeat = wx.TextCtrl(self, size=(182+25, 20), pos=(170, 215-2))
        self.textctrlRepeat.SetEditable(False)
        self.textctrlRepeat.SetForegroundColour((156, 156, 156)) #srsly?
        self.statictextSpacers = wx.StaticText(self, label="Spacers:", pos=(125, 240))
        self.listSpacers = wx.ListBox(self, id=wx.ID_ANY, size=(227+25, 90), pos=(125, 260))

        #genes
        self.staticboxGenes = wx.StaticBox(self, -1, "Genes", size=(379, 155), pos=(5, 365))
        self.listctrlGenes = wx.ListCtrl(self, style=wx.LC_REPORT, size=(363, 125), pos=(14, 385))
        self.listctrlGenes.InsertColumn(0, 'Subtype')
        self.listctrlGenes.InsertColumn(1, 'Name')
        self.listctrlGenes.InsertColumn(2, 'Symbol')
        self.listctrlGenes.InsertColumn(3, 'Start')
        self.listctrlGenes.InsertColumn(4, 'End')
        self.listctrlGenes.SetColumnWidth(0, 55)
        self.listctrlGenes.SetColumnWidth(1, 109)
        self.listctrlGenes.SetColumnWidth(2, 68)
        self.listctrlGenes.SetColumnWidth(3, 55)
        self.listctrlGenes.SetColumnWidth(4, 55)

        #bind keypress
        self.listAccessions.Bind(wx.EVT_KEY_DOWN, self.OnKeyPressOrganismAcc)
        self.listArrays.Bind(wx.EVT_KEY_DOWN, self.OnKeyPressArraysAcc)
        self.listSpacers.Bind(wx.EVT_KEY_DOWN, self.OnKeyPressSpacers)
        self.listctrlGenes.Bind(wx.EVT_KEY_DOWN, self.OnKeyPressGenes)
        
    def OnClickListGenomes(self, event):
        #get selected organism
        dataOrganismName = self.listGenome.Items[self.listGenome.GetSelection()] #WOW: Python needs this cached or finding the tuple takes ~5 seconds!
        dataOrganismTuple = [element for element in self.crisprIndex if element[2] == dataOrganismName][0]

        #update organism basic information
        self.textctrlName.SetValue(dataOrganismTuple[2])
        self.textctrlRefSeq.SetValue(dataOrganismTuple[0])
        self.textctrlTaxID.SetValue(dataOrganismTuple[1])
        self.listAccessions.Items = dataOrganismTuple[3].split(",")

        #update array data
        arrayData = self.dictCRISPRs[dataOrganismTuple[1]]
        accessions = arrayData.keys()
        accessions.sort()
        self.listArrays.Items = accessions
        self.textctrlRepeat.SetValue("")
        self.listSpacers.Items = []

        #update gene data
        self.listctrlGenes.ClearAll()
        self.listctrlGenes.InsertColumn(0, 'Subtype')
        self.listctrlGenes.InsertColumn(1, 'Name')
        self.listctrlGenes.InsertColumn(2, 'Symbol')
        self.listctrlGenes.InsertColumn(3, 'Start')
        self.listctrlGenes.InsertColumn(4, 'End')
        self.listctrlGenes.SetColumnWidth(0, 55)
        self.listctrlGenes.SetColumnWidth(1, 109)
        self.listctrlGenes.SetColumnWidth(2, 68)
        self.listctrlGenes.SetColumnWidth(3, 55)
        self.listctrlGenes.SetColumnWidth(4, 55)

        foundData = False
        for key in self.casGenes.keys():
            entry = self.casGenes[key]

            if dataOrganismName == entry[0]:
                geneDict = entry[1]
                subtypes = geneDict.keys()

                for type in subtypes:
                    for gene in geneDict[type].keys():
                        foundData = True
                        pos = self.listctrlGenes.InsertStringItem(0, type)
                        self.listctrlGenes.SetStringItem(pos, 1, gene)
                        self.listctrlGenes.SetStringItem(pos, 2, geneDict[type][gene][0])
                        self.listctrlGenes.SetStringItem(pos, 3, geneDict[type][gene][1][0][1])
                        self.listctrlGenes.SetStringItem(pos, 4, geneDict[type][gene][1][0][0])
                        
                break

        if not foundData:
            pos = self.listctrlGenes.InsertStringItem(0, "No")
            self.listctrlGenes.SetStringItem(pos, 1, "JCVI")
            self.listctrlGenes.SetStringItem(pos, 2, "Data")

    def OnClickArray(self, event):
        #get selected organism
        dataOrganismName = self.listGenome.Items[self.listGenome.GetSelection()] #WOW: Python needs this cached or finding the tuple takes ~5 seconds!
        dataOrganismTuple = [element for element in self.crisprIndex if element[2] == dataOrganismName][0]

        arrays = self.dictCRISPRs[dataOrganismTuple[1]]
        array = arrays[self.listArrays.Items[self.listArrays.GetSelection()]]

        self.textctrlRepeat.SetValue(array[0])
        self.listSpacers.Items = array[1]

    def OnKeyPressOrganismAcc(self, event):
        if event.GetKeyCode() == wx.WXK_CONTROL:
            clipdata = wx.TextDataObject()
            clipdata.SetText(self.listAccessions.Items[self.listAccessions.GetSelection()])
            wx.TheClipboard.Open()
            wx.TheClipboard.SetData(clipdata)
            wx.TheClipboard.Close()

    def OnKeyPressArraysAcc(self, event):
        if event.GetKeyCode() == wx.WXK_CONTROL:
            clipdata = wx.TextDataObject()
            clipdata.SetText(self.listArrays.Items[self.listArrays.GetSelection()])
            wx.TheClipboard.Open()
            wx.TheClipboard.SetData(clipdata)
            wx.TheClipboard.Close()

    def OnKeyPressSpacers(self, event):
        if event.GetKeyCode() == wx.WXK_CONTROL:
            clipdata = wx.TextDataObject()
            clipdata.SetText(self.listSpacers.Items[self.listSpacers.GetSelection()])
            wx.TheClipboard.Open()
            wx.TheClipboard.SetData(clipdata)
            wx.TheClipboard.Close()

    def OnKeyPressGenes(self, event):
        if event.GetKeyCode() == wx.WXK_CONTROL:
            index = self.listctrlGenes.GetFirstSelected()
            text = self.listctrlGenes.GetItem(index, 0).GetText() + "    "
            text += self.listctrlGenes.GetItem(index, 1).GetText() + "    "
            text += self.listctrlGenes.GetItem(index, 2).GetText() + "    "
            text += self.listctrlGenes.GetItem(index, 3).GetText() + "    "
            text += self.listctrlGenes.GetItem(index, 4).GetText()

            clipdata = wx.TextDataObject()
            clipdata.SetText(text)
            wx.TheClipboard.Open()
            wx.TheClipboard.SetData(clipdata)
            wx.TheClipboard.Close()

class PanelSpacers(wx.Panel):
    def BindGenomes(self):
        #bind outside GUI
        #parent.Bind(wx.EVT_LISTBOX, self.OnClickListGenomes, listGenomes)
        self.parent.listGenomes.Bind(wx.EVT_LISTBOX, self.OnClickListGenomes)
        
    def __init__(self, parent, dimensions, location, listGenomes, crisprIndex, dictCRISPRs):
        wx.Panel.__init__(self, parent, size=dimensions, pos=location)

        self.parent = parent
        self.BindGenomes()
        self.listGenome = listGenomes
        self.crisprIndex = crisprIndex
        self.dictCRISPRs = dictCRISPRs

        #repeats
        self.statictextboxRepeat = wx.StaticBox(self, -1, "Repeat", size=(379, 170), pos=(5, 10))
        
        self.listctrlRepeats = wx.ListCtrl(self, style=wx.LC_REPORT, size=(361, 140), pos=(14, 30))
        self.Bind(wx.EVT_LIST_ITEM_SELECTED, self.OnClickListRepeats, self.listctrlRepeats)
        self.listctrlRepeats.InsertColumn(0, 'Accession')
        self.listctrlRepeats.InsertColumn(1, 'Sequence')
        self.listctrlRepeats.InsertColumn(2, 'Spacers')
        self.listctrlRepeats.SetColumnWidth(0, 84)
        self.listctrlRepeats.SetColumnWidth(1, 221)
        self.listctrlRepeats.SetColumnWidth(2, 52)

        #spacer length
        self.statictextboxSpacer = wx.StaticBox(self, -1, "Spacer", size=(379, 212), pos=(5, 185))

        self.statictextSpacerLength = wx.StaticText(self, label="Length: ", pos=(14, 205))
        self.textctrlSpacerLength = wx.TextCtrl(self, size=(321, 20),pos=(55, 202))
        self.textctrlSpacerLength.SetEditable(False)
        self.textctrlSpacerLength.SetForegroundColour((156, 156, 156)) #srsly?

        self.radioSeq = wx.RadioButton(self, -1, 'Enter Sequence:', pos=(14, 240), style=wx.RB_GROUP)
        self.radioAccession = wx.RadioButton(self, -1, 'NCBI Accession:', pos=(14, 340))
        self.radioFASTA = wx.RadioButton(self, -1, 'Select FASTA:', pos=(14, 370)) #off by 5

        self.Bind(wx.EVT_RADIOBUTTON, self.OnSelectRadioFASTA, self.radioFASTA)

        self.textctrlSeq = wx.TextCtrl(self, size=(261, 90), pos=(115, 240), style=wx.TE_MULTILINE | wx.TE_DONTWRAP)
        self.textctrlAccession = wx.TextCtrl(self, size=(261, 20), pos=(115, 340))
        self.buttonSelectFASTA = wx.Button(self, label="Select File", pos=(115,365))
        self.Bind(wx.EVT_BUTTON, self.OnClickSelectFASTA, self.buttonSelectFASTA)
        self.fastaFilename = ""

        #toogle hairpinnning
        self.checkboxHairpinning = wx.CheckBox(self, label='Filter hairpinning (experimental)', pos=(14, 430))

        #compute button
        self.buttonComputeSpacers = wx.Button(self, label="Compute Spacers", pos=(142,450))
        self.Bind(wx.EVT_BUTTON, self.OnClickComputeSpacers, self.buttonComputeSpacers)

    def OnSelectRadioFASTA(self, event):
        self.fastaPath = ""

    def OnClickSelectFASTA(self, event):
        dialog = wx.FileDialog(self, "Select FASTA File", "", "", "FASTA files (*.fasta)|*.fasta|Text files (*.txt)|*.txt", wx.OPEN)

        if dialog.ShowModal() == wx.ID_OK:
            self.fastaPath = dialog.GetDirectory() + "\\" + dialog.GetFilename()

        dialog.Destroy()

    def OnClickListGenomes(self, event):
        #get selected organism
        dataOrganismName = self.listGenome.Items[self.listGenome.GetSelection()] #WOW: Python needs this cached or finding the tuple takes ~5 seconds!
        dataOrganismTuple = [element for element in self.crisprIndex if element[2] == dataOrganismName][0]
        dictOrganismLocus = self.dictCRISPRs[str(dataOrganismTuple[1])]
        
        #update repeats
        self.listctrlRepeats.DeleteAllItems()
        for i, key in enumerate(dictOrganismLocus.keys()):
            self.listctrlRepeats.InsertStringItem(i, key)
            self.listctrlRepeats.SetStringItem(i, 1, dictOrganismLocus[key][0])
            self.listctrlRepeats.SetStringItem(i, 2, str(len(dictOrganismLocus[key][1])))

        #clear spacer length
        self.textctrlSpacerLength.SetValue("")

    def OnClickListRepeats(self, event):
        #get selected organism
        dataOrganismName = self.listGenome.Items[self.listGenome.GetSelection()] #WOW: Python needs this cached or finding the tuple takes ~5 seconds!
        dataOrganismTuple = [element for element in self.crisprIndex if element[2] == dataOrganismName][0]
        dictOrganismLocus = self.dictCRISPRs[str(dataOrganismTuple[1])]

        self.textctrlSpacerLength.SetValue(str(len(dictOrganismLocus[self.listctrlRepeats.GetItem(self.listctrlRepeats.GetFirstSelected(), 0).GetText()][1][0])))

    def OnClickComputeSpacers(self, event):
        print "Computing spacers..."

        #determine organism
        if self.listGenome.GetSelection() == -1:
            print "    Please select an organism."
            return

        dataOrganismName = self.listGenome.Items[self.listGenome.GetSelection()] #WOW: Python needs this cached or finding the tuple takes ~5 seconds!
        dataOrganismTuple = [element for element in self.crisprIndex if element[2] == dataOrganismName][0]

        #determine repeat
        if self.listctrlRepeats.GetFirstSelected() == -1:
            print "    Please select a repeat."
            return
        
        dataRepeat = self.listctrlRepeats.GetItem(self.listctrlRepeats.GetFirstSelected(), 1).GetText()

        #determine spacer length.
        if self.textctrlSpacerLength.GetValue() == "":
            print "    Please select a repeat with known spacers."
            return
        dataSpacerLength = int(self.textctrlSpacerLength.GetValue())

        #determine sequence
        #TODO: enforce DNA sequence
        if self.radioSeq.GetValue():
            if self.textctrlSeq.GetValue() == "":
                print "    Please select a sequence. \n"
                return
            dataSequence = self.textctrlSeq.GetValue()
            dataSequence = "".join(dataSequence.split())
        elif self.radioAccession.GetValue():
            if self.textctrlAccession.GetValue() == "":
                print "    Please enter an accession number with version."
                return

            if not "." in self.textctrlAccession.GetValue():
                print "    Please list version with the accession number."
                return

            accession = self.textctrlAccession.GetValue()
            filename = accession + ".fasta"

            if not os.path.exists(filename):
                Entrez.email = ''
                fetch = Entrez.efetch(db='nucleotide', id=accession, rettype='fasta')
                file = open(filename, 'w')
                file.write(fetch.read())
                file.close()
                fetch.close()

            dataSequence = SeqIO.read(filename, "fasta").seq.tostring()
            
        else: #FASTA
            fileFASTA = open(self.fastaPath)
            dataSequence = SeqIO.read(fileFASTA, "fasta").seq.tostring()
            fileFASTA.close()

        #make sure the sequence is long enough
        if len(dataSequence) < dataSpacerLength:
            print "    Input sequence is too short. \n"
            return

        dataSequence = dataSequence.upper()

        #check if the input sequence might not be DNA
        for residue in dataSequence:
            #okay, so these are also single letter AAs, but close enough for now.
            if not (residue == "A" or residue == "T" or residue == "C" or residue == "G"):
                print "  Sequence is not nucleotides!"
                return

        print "    Organism: %s" % dataOrganismTuple[2]
        print "    Repeat: %s" % dataRepeat
        print "    Spacer Length: %s" % dataSpacerLength
        print "    Sequence: %s" % dataSequence

        #(sequence, repeat, organismUID, spacerLength)

            
        spacers = ArrayGeneration.getValidSpacers(Seq(dataSequence, IUPAC.unambiguous_dna),
                                                  Seq(dataRepeat, IUPAC.unambiguous_dna),
                                                  dataOrganismTuple[0],
                                                  dataSpacerLength,
                                                  self.checkboxHairpinning.GetValue())            


        dialog = DialogSpacersFinalize(self, spacers, dataRepeat)
        dialog.ShowModal()
        dialog.Destroy()
        
class MainFrame(wx.Frame):
    def __init__(self, parent, title):

        wx.Frame.__init__(self, parent, title=title, size=(800, 700), style=wx.DEFAULT_FRAME_STYLE ^ wx.RESIZE_BORDER)

        menubar = wx.MenuBar()
        self.SetMenuBar(menubar)

        menuView = wx.Menu()
        menubar.Append(menuView,"&View")
        entryInformation = menuView.Append(wx.ID_ANY, "Information")
        self.Bind(wx.EVT_MENU, self.OnViewInformation, entryInformation)
        entrySpacerCreator = menuView.Append(wx.ID_ANY,"&Spacer Creator")
        self.Bind(wx.EVT_MENU, self.OnViewSpacerCreator, entrySpacerCreator)

        menuHelp = wx.Menu()
        entryInformation = menuHelp.Append(wx.ID_ABOUT, "&About")
        self.Bind(wx.EVT_MENU, self.OnAboutAbout, entryInformation)
        menubar.Append(menuHelp, "&Help")

        menuSettings = wx.Menu()
        settingsInformation_RE = menuSettings.Append(wx.ID_ANY, "Spacer Generation: Restriction")
        self.Bind(wx.EVT_MENU, self.OnSpacerGenSettings_RE, settingsInformation_RE)
        settingsInformation_PAM = menuSettings.Append(wx.ID_ANY, "Spacer Generation: PAMs")
        self.Bind(wx.EVT_MENU, self.OnSpacerGenSettings_PAM, settingsInformation_PAM)
                
        menubar.Append(menuSettings, "&Settings")

        #genome panel
        self.panelGenomes = wx.Panel(self, size=(400, 529), pos=(0, 0))

        #load genome information
        self.crisprIndex = pickle.load(open("data/tax_info_filtered.txt"))
        self.dictCRISPRs = pickle.load(open("data/arrayDB.txt"))

        #list of available genomes
        self.statictextboxGenomes = wx.StaticBox(self.panelGenomes, -1, "Genomes", size=(385, 510), pos=(10, 10))
        self.statictextSearch = wx.StaticText(self.panelGenomes, label="Search:", pos=(20, 29))
        self.textctrlSearchGenomes = wx.TextCtrl(self.panelGenomes, size=(322, 20),pos=(64, 27))
        self.Bind(wx.EVT_TEXT, self.OnChangeSearchGenomes, self.textctrlSearchGenomes)

        self.listGenomes = wx.ListBox(self.panelGenomes, id=wx.ID_ANY, size=(366, 455), pos=(20, 55))
        self.OnChangeSearchGenomes(None)

        self.panelSpacers = PanelSpacers(self, dimensions=(795, 529), location=(400,0), listGenomes=self.listGenomes, crisprIndex=self.crisprIndex, dictCRISPRs=self.dictCRISPRs)
        self.panelSpacers.Hide()
        self.panelBrowser = PanelBrowser(self, dimensions=(795, 529), location=(400,0), listGenomes=self.listGenomes, crisprIndex=self.crisprIndex, dictCRISPRs=self.dictCRISPRs)
        
        self.OnViewSpacerCreator(None)
        
        #console
        self.panelConsole = wx.Panel(self, size=(795, 123), pos=(0, 529))
        self.panelConsole.statictextConsole = wx.StaticText(self.panelConsole, label="Console:", pos=(10, 0))
        self.panelConsole.textctrlConsole = wx.TextCtrl(self.panelConsole, size=(775, 100), pos=(10, 15), style=wx.TE_MULTILINE | wx.TE_DONTWRAP)
        self.panelConsole.textctrlConsole.SetFont(wx.Font(pointSize=10, family=wx.FONTFAMILY_MODERN, style=wx.NORMAL, weight=wx.FONTWEIGHT_NORMAL))
        self.panelConsole.textctrlConsole.SetEditable(False)
        self.panelConsole.textctrlConsole.SetForegroundColour((156, 156, 156)) #srsly?
        self.panelConsole.textctrlConsole.AppendText("Ready.\n")

        #redirect terminal text
        target = TextRedirector(self.panelConsole.textctrlConsole)
        sys.stdout = target
        
        self.Show(True)

    def OnAboutAbout(self,e):
        dialog = wx.MessageDialog(self, "CRISPRstudio: Spacer Selection and Sequence Repository", "About CRISPRstudio", wx.OK)
        dialog.ShowModal()
        dialog.Destroy()

    def OnSpacerGenSettings_RE(self, e):
        dialog = DialogSpacerGenSettings_RE(self)
        dialog.ShowModal()
        dialog.Destroy()
        
    def OnSpacerGenSettings_PAM(self, e):
        dialog = DialogSpacerGenSettings_PAM(self)
        dialog.ShowModal()
        dialog.Destroy()
        
    def OnChangeSearchGenomes(self, event):
        self.listGenomes.Items = [element[2] for element in self.crisprIndex if UtilGeneral.SAT(self.textctrlSearchGenomes.GetValue(), element[2])]

    def OnViewInformation(self, e):
        self.panelBrowser.Show()
        self.panelBrowser.BindGnomes()
        self.panelSpacers.Hide()

    def OnViewSpacerCreator(self,e):
        self.panelSpacers.Show()
        self.panelSpacers.BindGenomes()
        self.panelBrowser.Hide()

################################## FUNCTIONS ###################################

def testArrayGeneration():
    seq = Seq("CCATGGCCAACACTTGTCACTACTTTCGGT".lower(), IUPAC.unambiguous_dna)#DA
    repeat = Seq("TTTATCCCCGCTGGCGCGGGGAACTC", IUPAC.unambiguous_dna)
    organism = 57791
    spacerLength = 30

    print ArrayGeneration.getValidSpacers(seq, repeat, organism, spacerLength)

if __name__ == "__main__":


    try:
        last_update = open(Config.CRISPRDB_FILE_PREFIX + "last_update.txt").read()
        print "Last updated files on " + last_update
        update_choice = raw_input("Update files (y/n)? ")
        if update_choice == 'y':
            Updater.main()
    except IOError:
        print "Updater has not been run"
        update_choice = raw_input("Update files (y/n)? ")
        if update_choice == 'y':
            Updater.main()
        else:
            exit()

    app = wx.App(False)
    frame = MainFrame(None, "CRISPRstudio")
    app.MainLoop()

