from tkinter import *
from tkinter.ttk import *
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk) 
import os
import sys
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import mander
import steel

class Gui:

    def __init__(self):
        self.root = Tk()
        self.root.geometry("1300x900")
        self.root.title("CUMBIA")

        # Show the entries for circular section calculations
        self.button1 = Button(self.root, text = 'Run Analysis',command=self.run)
        self.button1.grid(row=0,column=0, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)
        self.button2 = Button(self.root, text = 'Quit',command=self._quit)
        self.button2.grid(row=0,column=1, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)

        # Set up the Combobox
        self.label1 = Label(self.root, text = "Memeber Section Type")
        self.label1.grid(row = 2, column = 0, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)
        self.box1 = Combobox(self.root)
        self.box1.grid(row = 2, column = 1, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)
        self.box1['values'] = ['Circular']
        self.box1.current(0)

        # The Entry to be shown
        self.show_circ = False; self.Circ = {}
        self.show_rect = False; self.Rect = {}
        
        # The dictionary containts data to plot figures
        self.figures = {}

        # Check the selection in 100 ms
        self.root.after(100, self.check_for_section)

    def _quit(self):
        self.root.quit()     # stops mainloop
        self.root.destroy()  # this is necessary on Windows to prevent
                        # Fatal Python Error: PyEval_RestoreThread: NULL tstate

    def run(self):
        section = self.box1.get()
        if section == 'Circular':
            e77 = Button(self.root, text = 'Plot',command=self.plot)
            e77.grid(row = 1,column=0, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)
            e78 = Combobox(self.root)
            e78.grid(row = 1, column = 1, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)
            e78['values'] = ['1: Stress-Strain Relation for Concrete', \
                             '2: Stress-Strain Relation for Reinforcing Steel',\
                             '3: Moment - Curvature Relation',\
                             '4: Moyer - Kowalsky Buckling Model',\
                             '5: Berry - Eberhard Buckling Model',\
                             '6: Force - Displacement Relation',\
                             '7: Interaction Diagram',\
                             '8: Concrete Strain - (M-N)',\
                             '9: Steel Strain - (M-N)']
            e78.current(5)   
            self.Circ['figures'] = e78
            self.Circ['plot'] = e77
            
            self.run_circ()
            self.plot()

    def check_for_section(self):
        '''Checks the section type, and shows relevant entries '''

        # Get the value of the Combobox
        section = self.box1.get()

        # If the value is equal to "Custom" and show_field is set to False
        if section == 'Circular' and not self.show_circ:

            # Set show_field to True and pack() the custom entry field
            self.show_circ = True
            self.show_rect = False
            
            self.circ_input()

        # If the value DOESNT equal "Custom"
        elif section == 'Rectangular' and not self.show_rect:
            self.show_circ = False
            self.show_rect = True   
            print(section)

            # Destroy the circ input
            for obj in self.Circ['gui_objects']:
                obj.destroy()

        # Call this method again to keep checking the selection box
        self.root.after(100, self.check_for_section)

    def plot(self):
        plt.close('all')
        try:
            # destroy all widgets from frame
            for widget in self.frame.winfo_children():
               widget.destroy()
        
            # this will clear frame and frame will be empty
            # if you want to hide the empty panel then
            self.frame.pack_forget()
        except:
            self.frame = Frame(self.root)
        self.frame.grid(row=1, column=3, columnspan=50, rowspan=50)
        fig_no = self.Circ['figures'].get()
        plt.ioff()
        fig, ax = plt.subplots(figsize=(8, 6), dpi = 100)
        if fig_no[0] == '1':
            ax.fill_between(self.figures['1']['ec'], 0, self.figures['1']['fc'], edgecolor = 'black', facecolor='green', interpolate=True, alpha = 0.5, label = 'Confined Concrete')
            ax.fill_between(self.figures['1']['ecun'], 0, self.figures['1']['fcun'], edgecolor = 'black', facecolor='blue', interpolate=True, alpha = 0.5, label = 'Unconfined Concrete')
            ax.set_ylabel('Stress [MPa]')
            ax.set_xlabel('Strain')
            ax.set_xlim([0,1.05*self.figures['1']['ec'][-3]])
            ax.set_ylim([0,1.05*np.max(self.figures['1']['fc'])])
            ax.legend(loc = 'upper right')
            ax.grid(True)
            ax.set_title('Stress-Strain Relation for Concrete')

        elif fig_no[0] == '2':
            ax.fill_between(self.figures['2']['es'], 0, self.figures['2']['fs'], edgecolor = 'black', facecolor='Red', interpolate=True, alpha = 0.5, label = 'Reinforcing Steel')
            ax.set_ylabel('Stress [MPa]')
            ax.set_xlabel('Strain')
            ax.set_xlim([1.05*self.figures['2']['es'][3],1.05*self.figures['2']['es'][-3]])
            ax.set_ylim([1.05*self.figures['2']['fs'][3],1.05*np.max(self.figures['2']['fs'])])
            ax.legend(loc = 'upper left')
            ax.grid(True)
            ax.set_title('Stress-Strain Relation for Reinforcing Steel')

        elif fig_no[0] == '3':
            ax.plot(self.figures['3']['curvbilin'],self.figures['3']['mombilin'],c='red')
            ax.plot(self.figures['3']['curv'],self.figures['3']['mom'],c='blue',linestyle='--')
            ax.grid(True)
            ax.set_xlabel('Curvature [m$^{-1}$]',fontsize=16)
            ax.set_ylabel('Moment [kN.m]',fontsize=16)
            ax.set_title('Moment - Curvature Relation',fontsize=16)

        elif fig_no[0] == '4':
            if self.figures['4']['opt'] == 1:
                ax.plot(self.figures['4']['failCuDuMK'], self.figures['4']['failss'],'m', marker = '.', MarkerEdgeColor = 'k', MarkerFaceColor = 'g', MarkerSize = 12, label = None)
                
            ax.plot(self.figures['4']['CuDu'],-self.figures['4']['steelstrain'], c = 'red', label = 'Column strain ductility behavior')
            ax.plot(self.figures['4']['CuDu'],self.figures['4']['esfl'],c = 'blue', linestyle='--', label = 'Flexural Tension Strain')
            ax.grid(True)
            ax.set_xlabel('Curvature Ductility',fontsize=16)
            ax.set_ylabel('Steel Tension Strain',fontsize=16)
            ax.legend()
            ax.set_title('Moyer - Kowalsky Buckling Model',fontsize=16)

        elif fig_no[0] == '5':
            if self.figures['5']['opt'] == 1:
                ax.plot(self.figures['5']['failCuDuBE'],self.figures['5']['failplrot'], linestyle='', marker = '.', MarkerEdgeColor = 'k', MarkerFaceColor = 'g', MarkerSize = 12, label = 'Buckling')
        
            ax.plot(self.figures['5']['CuDu'],self.figures['5']['rotb']*np.ones(len(self.figures['5']['CuDu'])), c = 'red', label = 'Plastic Rotation for Buckling')
            ax.plot(self.figures['5']['CuDu'],self.figures['5']['plrot'],c = 'blue', linestyle='--', label = 'Plastic Rotation')    
            ax.set_xlabel('Curvature Ductility',fontsize=16)
            ax.set_ylabel('Plastic Rotation',fontsize=16)
            ax.legend()
            ax.set_title('Berry - Eberhard Buckling Model',fontsize=16)
            ax.grid(True)

        elif fig_no[0] == '6':
            fig, ax = plt.subplots(figsize=(8, 7.5), dpi = 100)
            # Shink current axis by 20%
            box = ax.get_position()
            ax.set_position([box.x0, box.y0+box.height*0.2, box.width, box.height * 0.8])

            if self.figures['6']['opt'] == 1: 
                ax.plot(self.figures['6']['faildispl'],self.figures['6']['failforce'], c = 'magenta', linestyle = '', marker = '.', MarkerEdgeColor = 'black', MarkerFaceColor = 'red', markersize=15, label = 'shear failure')
                ax.plot(self.figures['6']['buckldispl'],self.figures['6']['bucklforce'], c = 'red', linestyle = '', marker = '*', MarkerEdgeColor = 'black', MarkerFaceColor = 'red', markersize=15, label = 'buckling (M & K)') 

            elif self.figures['6']['opt'] == 2: 
                ax.plot(self.figures['6']['faildispl'],self.figures['6']['failforce'], c = 'magenta', linestyle = '', marker = '.', MarkerEdgeColor = 'black', MarkerFaceColor = 'red', markersize=15, label = 'shear failure')
                ax.plot(self.figures['6']['buckldispl'],self.figures['6']['bucklforce'], c = 'red', linestyle = '', marker = '*', MarkerEdgeColor = 'black', MarkerFaceColor = 'red', markersize=15, label = 'buckling (M & K)') 
                ax.plot(self.figures['6']['buckldisplBE'],self.figures['6']['bucklforceBE'], c = 'red', linestyle = '', marker = 's', MarkerEdgeColor = 'black', MarkerFaceColor = 'red', markersize=15, label = 'buckling (B & E)')

            elif self.figures['6']['opt'] == 3: 
                ax.plot(self.figures['6']['faildispl'],self.figures['6']['failforce'], c = 'magenta', linestyle = '', marker = '.', MarkerEdgeColor = 'black', MarkerFaceColor = 'red', markersize=15, label = 'shear failure')
                ax.plot(self.figures['6']['buckldisplBE'],self.figures['6']['bucklforceBE'], c = 'red', linestyle = '', marker = 's', MarkerEdgeColor = 'black', MarkerFaceColor = 'red', markersize=15, label = 'buckling (B & E)')
            
            elif self.figures['6']['opt'] == 4: 
                ax.plot(self.figures['6']['faildispl'],self.figures['6']['failforce'], c = 'magenta', linestyle = '', marker = '.', MarkerEdgeColor = 'black', MarkerFaceColor = 'red', markersize=15, label = 'shear failure')

            elif self.figures['6']['opt'] == 5: 
                ax.plot(self.figures['6']['buckldispl'],self.figures['6']['bucklforce'], c = 'red', linestyle = '', marker = '*', MarkerEdgeColor = 'black', MarkerFaceColor = 'red', markersize=15, label = 'buckling (M & K)') 

            elif self.figures['6']['opt'] == 6: 
                ax.plot(self.figures['6']['buckldispl'],self.figures['6']['bucklforce'], c = 'red', linestyle = '', marker = '*', MarkerEdgeColor = 'black', MarkerFaceColor = 'red', markersize=15, label = 'buckling (M & K)') 
                ax.plot(self.figures['6']['buckldisplBE'],self.figures['6']['bucklforceBE'], c = 'red', linestyle = '', marker = 's', MarkerEdgeColor = 'black', MarkerFaceColor = 'red', markersize=15, label = 'buckling (B & E)')
                
            elif self.figures['6']['opt'] == 7: 
                ax.plot(self.figures['6']['buckldisplBE'],self.figures['6']['bucklforceBE'], c = 'red', linestyle = '', marker = 's', MarkerEdgeColor = 'black', MarkerFaceColor = 'red', markersize=15, label = 'buckling (B & E)')
        
            ax.plot(self.figures['6']['displ'],self.figures['6']['Force'], c = 'black', linestyle = '-', label = 'total response')
            ax.plot(self.figures['6']['displbilin'],self.figures['6']['forcebilin'], c = 'blue', linestyle = '--', label = 'bilinear approximation')
            ax.plot(self.figures['6']['displ'],self.figures['6']['V'], c = 'red', linestyle = ':', label = 'shear capacity (assessment)')
            ax.plot(self.figures['6']['displ'],self.figures['6']['Vd'], c = 'magenta', linestyle = ':', label = 'shear capacity (design)')
            ax.fill_between(self.figures['6']['displ'],0,self.figures['6']['Force'], edgecolor = 'black', facecolor=(0.4, 0.8, 0.6), interpolate=True, alpha = 0.4, label = 'ultimate zone')
            ax.fill_between(self.figures['6']['displ'][self.figures['6']['pointsdam']],0,self.figures['6']['Force'][self.figures['6']['pointsdam']], edgecolor = 'black', facecolor=(0.4, 0.6, 0.6), interpolate=True, alpha = 0.4, label = 'damage control zone')
            ax.fill_between(self.figures['6']['displ'][self.figures['6']['pointsser']],0,self.figures['6']['Force'][self.figures['6']['pointsser']], edgecolor = 'black', facecolor=(0.4, 0.4, 0.6), interpolate=True, alpha = 0.4, label = 'serviceability zone')
            ax.grid(True)
            ax.set_xlabel('Displacement [m]', fontsize=16)
            ax.set_ylabel('Force [kN]', fontsize=16)
            ax.set_title('Force - Displacement Relation', fontsize=16)
            ax.legend(bbox_to_anchor=(-0.1, -0.3), loc='lower left', frameon=False, ncol = 3)
            ylims = ax.get_ylim(); ax.set_ylim((0,ylims[1]))
            xlims = ax.get_xlim(); ax.set_xlim((0,xlims[1]))

        elif fig_no[0] == '7':
            ax.plot(self.figures['7']['Mni'],self.figures['7']['PPn']/1000, c = 'red', linestyle = '-', marker = 'o', MarkerEdgeColor = 'black', MarkerFaceColor = 'red', markersize=6, label = 'Interaction Diagram')
            ax.plot(self.figures['7']['MnL'],self.figures['7']['PPL']/1000, c = 'blue', linestyle = '--', marker = 's', MarkerEdgeColor = 'black', MarkerFaceColor = 'blue', markersize=6, label = 'Approximation for NLTHA')
            ax.set_xlabel('Moment [kN.m]', fontsize = 16)
            ax.set_ylabel('Axial Load [kN]', fontsize = 16)
            ax.set_title('Interaction Diagram', fontsize=16)
            ax.legend()
            ax.grid(True)

        elif fig_no[0] == '8':
            ax.plot(self.figures['8']['PPn'][1:-1]/1000, self.figures['8']['eci'], c = 'blue', linestyle = '--', marker = 's', MarkerEdgeColor = 'black', MarkerFaceColor = 'blue', markersize=6)
            ax.set_xlabel('Axial Force [kN]', fontsize = 16)
            ax.set_ylabel('Concrete Strain', fontsize = 16)
            ax.grid(True)

        elif fig_no[0] == '9':
            plt.plot(self.figures['9']['PPn'][1:-1]/1000, self.figures['9']['esi'], c = 'blue', linestyle = '--', marker = 's', MarkerEdgeColor = 'black', MarkerFaceColor = 'blue', markersize=6)
            ax.set_xlabel('Axial Force [kN]', fontsize = 16)
            ax.set_ylabel('Steel Strain', fontsize = 16)
            ax.grid(True)

        # creating the Tkinter canvas 
        # containing the Matplotlib figure 
        canvas = FigureCanvasTkAgg(fig, 
                                   master = self.frame)   
        canvas.draw() 
      
        # placing the canvas on the Tkinter window 
        get_widz = canvas.get_tk_widget()
        get_widz.pack()
      
        # creating the Matplotlib toolbar 
        toolbar = NavigationToolbar2Tk(canvas, 
                                        self.frame) 
        toolbar.update() 
      
        # placing the toolbar on the Tkinter window 
        canvas.get_tk_widget().pack()

    def circ_input(self):
        # Diameter
        e1=Label(self.root,text="Diameter [mm]")
        e1.grid(row = 3, column = 0, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)
        self.Circ['D'] = StringVar(value="979")
        e2=Entry(self.root,textvariable = self.Circ['D'])
        e2.grid(row = 3, column = 1, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)

        # Cover to longitudinal bars 
        e3=Label(self.root,text ="Cover to longitudinal bars [mm]")
        e3.grid(row = 4, column = 0, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)
        self.Circ['clb'] = StringVar(value = "28") 
        e4=Entry(self.root, textvariable = self.Circ['clb']) 
        e4.grid(row = 4, column = 1, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)
        
        # Member length
        e5=Label(self.root,text="Member clear length [mm]")
        e5.grid(row = 5, column = 0, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)
        self.Circ['L']=StringVar(value = "1651")
        e6=Entry(self.root,textvariable = self.Circ['L'])
        e6.grid(row = 5, column = 1, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)           
        
        # Bending type
        e7=Label(self.root,text="Bending type")
        e7.grid(row = 6, column = 0, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)
        e8 = Combobox(self.root)
        e8.grid(row = 6, column = 1, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)
        e8['values'] = ['single', 'double']
        e8.current(0)   
        self.Circ['bending'] = e8

        # Ductility mode
        e9=Label(self.root,text="Ductility mode")
        e9.grid(row = 7, column = 0, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)
        e10 = Combobox(self.root)
        e10.grid(row = 7, column = 1, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)
        e10['values'] = ['uniaxial', 'biaxial']
        e10.current(0)   
        self.Circ['ductilitymode'] = e10
 
        # Number of longitudinal bars
        e11=Label(self.root,text="Number of longitudinal bars")
        e11.grid(row = 8, column = 0, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)
        self.Circ['nbl']=StringVar(value = "30")
        e12=Entry(self.root,textvariable = self.Circ['nbl'])
        e12.grid(row = 8, column = 1, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)              

        # Long. bar diameter 
        e13=Label(self.root,text="Long. bar diameter [mm]")
        e13.grid(row = 9, column = 0, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)
        self.Circ['Dbl']=StringVar(value = "32")
        e14=Entry(self.root,textvariable = self.Circ['Dbl'])
        e14.grid(row = 9, column = 1, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)   

        # Diameter of transverse reinf.
        e15=Label(self.root,text="Diameter of transverse reinf. [mm]")
        e15.grid(row = 10, column = 0, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)
        self.Circ['Dh']=StringVar(value = "16")
        e16=Entry(self.root,textvariable = self.Circ['Dh'])
        e16.grid(row = 10, column = 1, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)   

        # Type of transverse reinf.
        e17=Label(self.root,text="Type of transverse reinf.")
        e17.grid(row = 11, column = 0, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)
        e18 = Combobox(self.root)
        e18.grid(row = 11, column = 1, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)
        e18['values'] = ['hoops', 'spirals']
        e18.current(1)   
        self.Circ['stype'] = e18

        # Spacing of transverse steel
        e19=Label(self.root,text="Spacing of transverse steel [mm]")
        e19.grid(row = 12, column = 0, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)
        self.Circ['s']=StringVar(value = "100")
        e20=Entry(self.root,textvariable = self.Circ['s'])
        e20.grid(row = 12, column = 1, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)
        
        # Axial load
        e21=Label(self.root,text="Axial load, -tens, +comp [kN]")
        e21.grid(row = 13, column = 0, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)
        self.Circ['P']=StringVar(value = "1325")
        e22=Entry(self.root,textvariable = self.Circ['P'])
        e22.grid(row = 13, column = 1, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)  

        # Concrete model
        e23=Label(self.root,text="Concrete type")
        e23.grid(row = 14, column = 0, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)
        e24 = Combobox(self.root)
        e24.grid(row = 14, column = 1, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)
        e24['values'] = ['Normal weight', 'Light weight']
        e24.current(0)
        self.Circ['Concrete'] = e24                

        # Steel model
        e25=Label(self.root,text="Reinforcing steel model")
        e25.grid(row = 15, column = 0, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)
        e26 = Combobox(self.root)
        e26.grid(row = 15, column = 1, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)
        e26['values'] = ['King', 'Raynor']
        e26.current(0)
        self.Circ['Steel'] = e26 

        # Concrete compressive strength
        e27=Label(self.root,text="Concrete compressive strength [MPa]")
        e27.grid(row = 16, column = 0, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)
        self.Circ['fpc']=StringVar(value = "{:.2f}".format(1.3*28))
        e28=Entry(self.root,textvariable = self.Circ['fpc'])
        e28.grid(row = 16, column = 1, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2) 
        
        # Concrete modulus of elasticity
        e29=Label(self.root,text="Concrete modulus of elasticity [MPa]\nInsert 0 for Mander model")
        e29.grid(row = 17, column = 0, columnspan = 1, rowspan = 2, sticky=W, pady = 2, padx = 2)
        self.Circ['Ec'] = StringVar(value = '0')
        e30=Entry(self.root,textvariable = self.Circ['Ec'])
        e30.grid(row = 17, column = 1, columnspan = 1, rowspan = 2, sticky=W, pady = 2, padx = 2)

        # Unconfined concrete strain at peak stress (usually 0.002 for normal weight or 0.004 for lightweight)
        e31=Label(self.root,text="Unconfined concrete \nStrain at peak stress")
        e31.grid(row = 19, column = 0, columnspan = 1, rowspan = 2, sticky=W, pady = 2, padx = 2)
        self.Circ['eco'] = StringVar(value = '0.002')
        e32=Entry(self.root,textvariable = self.Circ['eco'])
        e32.grid(row = 19, column = 1, columnspan = 1, rowspan = 2, sticky=W, pady = 2, padx = 2)
        
        # Max transv. steel strain (usually ~0.10-0.15)
        e33=Label(self.root,text="Max transverse steel strain")
        e33.grid(row = 21, column = 0, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)
        self.Circ['esm'] = StringVar(value = '0.10')
        e34=Entry(self.root,textvariable = self.Circ['esm'])
        e34.grid(row = 21, column = 1, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)

        # Max uncon. conc. strain (usually 0.0064)
        e35=Label(self.root,text="Max unconfined concrete strain")
        e35.grid(row = 22, column = 0, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)
        self.Circ['espall'] = StringVar(value = '0.0064')
        e36=Entry(self.root,textvariable = self.Circ['espall'])
        e36.grid(row = 22, column = 1, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)

        # Long steel yielding stress (MPa)
        e37=Label(self.root,text="Long steel yielding stress [MPa]")
        e37.grid(row = 23, column = 0, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)
        self.Circ['fy'] = StringVar(value = "{:.2f}".format(1.1*414))
        e38=Entry(self.root,textvariable = self.Circ['fy'])
        e38.grid(row = 23, column = 1, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)

        # Transverse steel yielding stress (MPa)
        e39=Label(self.root,text="Transverse steel yielding stress [MPa]")
        e39.grid(row = 24, column = 0, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)
        self.Circ['fyh'] = StringVar(value = "{:.2f}".format(1.1*414))
        e40=Entry(self.root,textvariable = self.Circ['fyh'])
        e40.grid(row = 24, column = 1, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)

        # Steel modulus of elasticity
        e41=Label(self.root,text="Steel modulus of elasticity [MPa]")
        e41.grid(row = 25, column = 0, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)
        self.Circ['Es'] = StringVar(value = '200000')
        e42=Entry(self.root,textvariable = self.Circ['Es'])
        e42.grid(row = 25, column = 1, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)

        # Long steel max stress (MPa)
        e43=Label(self.root,text="Longitudinal steel max stress [MPa]")
        e43.grid(row = 26, column = 0, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)
        self.Circ['fsu'] = StringVar(value = "{:.2f}".format(1.4*float(self.Circ['fy'].get())))
        e44=Entry(self.root,textvariable = self.Circ['fsu'])
        e44.grid(row = 26, column = 1, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)

        # Long steel strain for strain hardening (usually 0.008)
        e45=Label(self.root,text="Longitudinal steel strain for strain hardening")
        e45.grid(row = 27, column = 0, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)
        self.Circ['esh'] = StringVar(value = '0.008')
        e46=Entry(self.root,textvariable = self.Circ['esh'])
        e46.grid(row = 27, column = 1, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)
        
        # Longtudinal steel maximum strain (usually ~0.10-0.15)
        e47=Label(self.root,text="Longitudinal steel strain for strain hardening")
        e47.grid(row = 28, column = 0, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)
        self.Circ['esu'] = StringVar(value = '0.10')
        e48=Entry(self.root,textvariable = self.Circ['esu'])
        e48.grid(row = 28, column = 1, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)

        # Slope of the yield plateau (MPa)
        e49=Label(self.root,text="Slope of the yield plateau [MPa]")
        e49.grid(row = 29, column = 0, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)
        self.Circ['Ey'] = StringVar(value = '350')
        e50=Entry(self.root,textvariable = self.Circ['Ey'])
        e50.grid(row = 29, column = 1, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)

        # Slope of the yield plateau (MPa)
        e51=Label(self.root,text="Parameter defining strain hardening\ncurve in the Raynor model")
        e51.grid(row = 30, column = 0, columnspan = 1, rowspan = 2, sticky=W, pady = 2, padx = 2)
        self.Circ['C1'] = StringVar(value = '3.5')
        e52=Entry(self.root,textvariable = self.Circ['C1'])
        e52.grid(row = 30, column = 1, columnspan = 1, rowspan = 2, sticky=W, pady = 2, padx = 2)
        
        # strain limits for yield surface (interaction diagram);
        e53=Label(self.root,text="Concrete strain limit for yield surface\ninteraction diagram")
        e53.grid(row = 32, column = 0, columnspan = 1, rowspan = 2, sticky=W, pady = 2, padx = 2)
        self.Circ['csid'] = StringVar(value = '0.004')
        e54=Entry(self.root,textvariable = self.Circ['csid'])
        e54.grid(row = 32, column = 1, columnspan = 1, rowspan = 2, sticky=W, pady = 2, padx = 2)
        e55=Label(self.root,text="Steel strain limit for yield surface\ninteraction diagram")
        e55.grid(row = 34, column = 0, columnspan = 1, rowspan = 2, sticky=W, pady = 2, padx = 2)
        self.Circ['ssid'] = StringVar(value = '0.015')
        e56=Entry(self.root,textvariable = self.Circ['ssid'])
        e56.grid(row = 34, column = 1, columnspan = 1, rowspan = 2, sticky=W, pady = 2, padx = 2)

        # Deformation Limit States:
        e57=Label(self.root,text="Concrete serviceability strain")
        e57.grid(row = 36, column = 0, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)
        self.Circ['ecser'] = StringVar(value = '0.004')
        e58=Entry(self.root,textvariable = self.Circ['ecser'])
        e58.grid(row = 36, column = 1, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)
        e59=Label(self.root,text="Steel serviceability strain")
        e59.grid(row = 37, column = 0, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)
        self.Circ['esser'] = StringVar(value = '0.015')
        e60=Entry(self.root,textvariable = self.Circ['esser'])
        e60.grid(row = 37, column = 1, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)

        e61=Label(self.root,text="Concrete damage control strain")
        e61.grid(row = 38, column = 0, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)
        self.Circ['ecdam'] = StringVar(value = '0.018')
        e62=Entry(self.root,textvariable = self.Circ['ecdam'])
        e62.grid(row = 38, column = 1, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)
        e63=Label(self.root,text="Steel damage control strain")
        e63.grid(row = 39, column = 0, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)
        self.Circ['esdam'] = StringVar(value = '0.060')
        e64=Entry(self.root,textvariable = self.Circ['esdam'])
        e64.grid(row = 39, column = 1, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)      
        
        # temperature information (in case of freezing conditions)
        e65=Label(self.root,text="Temperature of the specimen in celsius")
        e65.grid(row = 40, column = 0, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)
        self.Circ['temp'] = StringVar(value = '40')
        e66=Entry(self.root,textvariable = self.Circ['temp'])
        e66.grid(row = 40, column = 1, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)
        
        # Constant to calculate strain penetration length
        e67=Label(self.root,text="Constant to calculate strain\npenetrationlengthusually 0.022 at\nambient temp. or 0.011 at -40C")
        e67.grid(row = 41, column = 0, columnspan = 1, rowspan = 3, sticky=W, pady = 2, padx = 2)
        self.Circ['kLsp'] = StringVar(value = '0.022')
        e68=Entry(self.root,textvariable = self.Circ['kLsp'])
        e68.grid(row = 41, column = 1, columnspan = 1, rowspan = 3, sticky=W, pady = 2, padx = 2)            
        
        # analysis control parameters
        e69=Label(self.root,text="Max number of iterations")
        e69.grid(row = 44, column = 0, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)
        self.Circ['itermax'] = StringVar(value = '1000')
        e70=Entry(self.root,textvariable = self.Circ['itermax'])
        e70.grid(row = 44, column = 1, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)
        e71=Label(self.root,text="Number of concrete layers")
        e71.grid(row = 45, column = 0, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)
        self.Circ['ncl'] = StringVar(value = '40')
        e72=Entry(self.root,textvariable = self.Circ['ncl'])
        e72.grid(row = 45, column = 1, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)
        e73=Label(self.root,text="Tolerance (x fpc x Ag)")
        e73.grid(row = 46, column = 0, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)
        self.Circ['tolerance'] = StringVar(value = '0.001')
        e74=Entry(self.root,textvariable = self.Circ['tolerance'])
        e74.grid(row = 46, column = 1, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)
        e75=Label(self.root,text="Strain step for default material models")
        e75.grid(row = 47, column = 0, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)
        self.Circ['dels'] = StringVar(value = '0.0001')
        e76=Entry(self.root,textvariable = self.Circ['dels'])
        e76.grid(row = 47, column = 1, columnspan = 1, rowspan = 1, sticky=W, pady = 2, padx = 2)

        self.Circ['gui_objects'] = []
        for i in range(1,77):
            self.Circ['gui_objects'].append(eval('e'+str(i)))
            
            
    def run_circ(self):
        
        name = 'Outputs'      # identifies actual work, the output file will be name.xls
        
        # section properties:
        D       = float(self.Circ['D'].get())        # section diameter (mm)
        clb     = float(self.Circ['clb'].get())            # cover to longitudinal bars (mm)
        
        # member properties
        
        L             = float(self.Circ['L'].get())       # member clear length (mm)
        bending       = self.Circ['bending'].get()    # single or double
        ductilitymode = self.Circ['ductilitymode'].get()  # biaxial or uniaxial
        
        # reinforcement details:
        
        nbl     = int(self.Circ['nbl'].get())             # number of longitudinal bars
        Dbl     = int(self.Circ['Dbl'].get())              # long. bar diameter (mm)   
        Dh      = int(self.Circ['Dh'].get())               # diameter of transverse reinf. (mm)
        stype   = self.Circ['stype'].get()     # 'spirals' or 'hoops'*
        s       = int(self.Circ['s'].get())            # spacing of transverse steel (mm)*
        
        # aplieed loads:
        
        P      =  float(self.Circ['P'].get())            # axial load kN (-) tension (+)compression
        
        # material models (input the 'name' of the file with the stress-strain relationship
        # to use the default models: Mander model for confined or unconfined  concrete type 'mc' or 'mu'.
        # For lightweight confined concrete type 'mclw' 
        # King model for the steel 'ks', Raynor model for steel 'ra':
            
        if self.Circ['Concrete'].get() == 'Normal weight':
            confined   = 'mc'
            unconfined = 'mu'
        elif self.Circ['Concrete'].get() == 'Light weight':
            confined   = 'mclw'
            unconfined = 'mulw'        

        if self.Circ['Steel'].get() == 'King':
            rebar      = 'ks'
        elif self.Circ['Steel'].get() == 'Raynor':
            rebar      = 'ra'
            
        # material properties 
        
        fpc     = float(self.Circ['fpc'].get())  # concrete compressive strength (MPa)
        Ec      = float(self.Circ['Ec'].get())        # concrete modulus of elasticity (MPa) or
                          # input 0 for automatic calculation using
                          # 5000(fpc)^0.5
        eco     = float(self.Circ['eco'].get())   # unconfined strain (usually 0.002 for normal weight or 0.004 for lightweight)*
        esm     = float(self.Circ['esm'].get())    # max transv. steel strain (usually ~0.10-0.15)*
        espall  = float(self.Circ['espall'].get())  # max uncon. conc. strain (usually 0.0064)
        
        fy      = float(self.Circ['fy'].get()) # long steel yielding stress (MPa)
        fyh     = float(self.Circ['fyh'].get()) # transverse steel yielding stress (MPa)
        Es      = float(self.Circ['Es'].get())  # steel modulus of elasticity
        fsu     = float(self.Circ['fsu'].get())  # long steel max stress (MPa)*
        esh     = float(self.Circ['esh'].get())  # long steel strain for strain hardening (usually 0.008)*
        esu     = float(self.Circ['esu'].get())    # long. steel maximum strain (usually ~0.10-0.15)*
        
        Ey     =  float(self.Circ['Ey'].get())     # slope of the yield plateau (MPa)
        C1     =  float(self.Circ['C1'].get())     # defines strain hardening curve in the Raynor model [2-6]
        
        # this information is used only if the default material models are selected
        
        # strain limits for yield surface (interaction diagram);
        
        csid = float(self.Circ['csid'].get())  # concrete
        ssid = float(self.Circ['ssid'].get())  # steel
        
        # Deformation Limit States:
        
        ecser = float(self.Circ['ecser'].get()); esser = float(self.Circ['esser'].get())   # concrete (ecser) and steel (esser) serviceability strain
        ecdam = float(self.Circ['ecdam'].get()); esdam = float(self.Circ['esdam'].get())   # concrete (ecser) and steel (esser) damage control strain
                                       # (to use the 2/3 of the ultimate concrete strain just tipe 'twth'
        # temperature information (in case of freezing conditions)
        temp  = float(self.Circ['temp'].get())           # temperature of the specimen in celsius
        kLsp  = float(self.Circ['kLsp'].get())        # constant to calculate Lsp = kLsp*fy*Dbl
                             # (usually 0.022 at ambient temp. or 0.011 at -40C)
        
        # analysis control parameters:
        itermax    = int(self.Circ['itermax'].get())       # max number of iterations
        ncl        = int(self.Circ['ncl'].get())         # # of concrete layers
        tolerance  = float(self.Circ['tolerance'].get())      # x fpc x Ag
        dels       = float(self.Circ['dels'].get())     # delta strain for default material models
        
        ### INPUT ENDS ###
        
        models_path = 'models' # user specified models

        # concrete modulus of elaticity
        if Ec == 0:
            Ec = 5000*(fpc**0.5)                  
        
        # tensile strength, considering temperature effect
        if temp < 0:
            fct = (1-0.0105*temp)*0.56*(fpc**0.5)
        elif temp >= 0:
            fct = 0.56*(fpc**0.5)
        
        eccr = fct/Ec                   # concrete strain for cracking
        
        Dsp = D - 2 * clb + Dh          # core diameter
        dcore = clb - Dh * 0.5          # distance to the core
        P      = P * 1000               # axial load in Newtons
        Ast    = nbl * 0.25 * np.pi * (Dbl**2) # total long. steel area mm2
        
        tcl = D / ncl                   # thickness of concrete layers
        yl  = tcl * np.arange(1,ncl+1)    # border distance conc. layer
        
        esser = -esser
        esdam = -esdam
        
        ecun, fcun = mander.circ_un(Ec,fpc,eco,espall,dels)
        ec, fc = mander.circ_conf(Ec,Ast,Dh,clb,s,fpc,fy,eco,esm,D,dels,stype)
        
        if unconfined == 'mu':
            ecun, fcun = mander.circ_un(Ec,fpc,eco,espall,dels)
        elif unconfined == 'mc':
            ecun, fcun = mander.circ_conf(Ec,Ast,Dh,clb,s,fpc,fy,eco,esm,D,dels,stype)
        elif unconfined == 'mclw':
            ecun, fcun = mander.circ_conflw(Ec,Ast,Dh,clb,s,fpc,fy,eco,esm,D,dels,stype)
        else:
            AUX = np.loadxtxt(os.path.join(models_path,unconfined+'txt'))
            ecun = AUX[:,0]
            fcun = AUX[:,1]
            
        if confined == 'mu':
            ec, fc = mander.circ_un(Ec,fpc,eco,espall,dels)
        elif confined == 'mc':
            ec, fc = mander.circ_conf(Ec,Ast,Dh,clb,s,fpc,fy,eco,esm,D,dels,stype)
        elif confined == 'mclw':
            ec, fc = mander.circ_conflw(Ec,Ast,Dh,clb,s,fpc,fy,eco,esm,D,dels,stype)
        else:
            AUX = np.loadxtxt(os.path.join(models_path,confined+'txt'))
            ec = AUX[:,0]
            fc = AUX[:,1]    
        
        if rebar == 'ks':
            es, fs = steel.king(Es,fy,fsu,esh,esu,dels)
        elif rebar == 'ra':
            es, fs = steel.raynor(Es,fy,fsu,esh,esu,dels,C1,Ey)
        else:
            AUX = np.loadxtxt(os.path.join(models_path,rebar+'txt'))
            es = AUX[:,0]
            fs = AUX[:,1]    
          
        ecu = ec[-1]                         # maximum strain confined concrete
        ecumander = ecu/1.5                  # ultimate strain predicted by the original mander model
        
        if ecdam == 'twth':
            ecdam = ecumander
        
        # vector with strains of confined concrete
        ec = np.append(-1e10, ec)
        ec = np.append(ec,ec[-1]+dels)
        ec = np.append(ec,1e10)
        # vector with stresses of confined concrete  
        fc = np.append(0, fc)
        fc = np.append(fc, 0)
        fc = np.append(fc, 0)
        
        # vector with strains of unconfined concrete
        ecun = np.append(-1e10, ecun)
        ecun = np.append(ecun,ecun[-1]+dels)
        ecun = np.append(ecun,1e10)
        # vector with stresses of unconfined concrete
        fcun = np.append(0, fcun)
        fcun = np.append(fcun, 0)
        fcun = np.append(fcun, 0)
        
        # maximum strain steel
        esu = es[-1]
        # vector with strains of the steel
        es = np.append(es, es[-1]+dels)
        es = np.append(es, 1e10)
        # vector with stresses of the steel
        fs = np.append(fs, 0)
        fs = np.append(fs, 0)
        
        esaux = np.zeros(len(es))
        fsaux = 0*esaux
        
        for i in range(len(es)):
            esaux[i] = es[len(es)-i-1]
            fsaux[i] = fs[len(fs)-i-1]
        
        # vector with strains of the steel
        es = np.append(-esaux, es[3:])
        # vector with stresses of the steel
        fs = np.append(-fsaux, fs[3:])
        
        self.figures['1'] = {'ec': ec, 'fc': fc, 'ecun':ecun, 'fcun':fcun}
        self.figures['2'] = {'es': es, 'fs': fs}
        
        # ============================== CONCRETE LAYERS ============================
        
        # add layers to consider unconfined concrete
        yl = np.append(np.append(yl, dcore), D-dcore)
        yl.sort()
        # confined concrete layers
        yc = yl - dcore;
        yc = yc[np.where((0 < yc)*(yc <= Dsp))[0]] 
        
        # total area of each layer
        Atemp = ((D/2)**2)*np.arccos(1-2*yl/D)-(D/2-yl)*((D*yl-yl**2)**0.5)
        Atc = Atemp-np.append(0,Atemp[:-1])
        
        # total area of each conf. layer
        Atemp = ((Dsp/2)**2)*np.arccos(1-2*yc/Dsp)-(Dsp/2-yc)*((Dsp*yc-yc**2)**0.5)
        Atcc = Atemp-np.append(0,Atemp[:-1])
        
        conclay = []
        k = 0
        for i in range(len(yl)):
            if yl[i]<=dcore or yl[i]>D-dcore:
                conclay.append([Atc[i], 0])
                
            if yl[i]>dcore and yl[i]<=D-dcore:
                conclay.append([Atc[i]-Atcc[k], Atcc[k]])
                k += 1
        
        conclay = np.asarray(conclay)
        yl.shape = (len(yl),1)
        ycenter = np.append(yl[0]/2, 0.5*(yl[:-1]+yl[1:])); ycenter.shape = (len(ycenter),1)
        
        # [center_layer|A_uncon|A_conf|d_top_layer]
        conclay = np.concatenate((ycenter,conclay,yl),axis = 1)
        
        # ================================    REBARS     =====================================
        
        Asb   = 0.25*np.pi*(Dbl**2)
        r     = 0.5*(D-2*clb-Dbl)
        theta = (2*np.pi/nbl)*np.arange(0,nbl)
        distld = (0.5*(D-2*r)+r*np.sin(theta)*np.tan(0.5*theta))
        distld.sort()         # y coordinate of each bar
        
        # =============================== CORRECTED AREAS ======================================
        
        # Substract the steel area
        for i in range(nbl):
            aux = np.where(yl > distld[i])[0][0]
            conclay[aux,2] = conclay[aux,2] - Asb
            if conclay[aux,2] < 0:
                print('decrease # of layers')
                sys.exit()
                
        # ============  Define vector (def) with the deformations in the top concrete ==================
        
        df = np.arange(1e-4,20*ecu,1e-4)
            
        if ecu > 0.0018:
            df = np.append(df[df<=16e-4],np.arange(18e-4,20*ecu,2e-4))
            
        if ecu > 0.0025:
            df = np.append(df[df<=20e-4],np.arange(25e-4,20*ecu,5e-4))
            
        if ecu > 0.006:
            df = np.append(df[df<=50e-4],np.arange(60e-4,20*ecu,10e-4))
            
        if ecu > 0.012: 
            df = np.append(df[df<=100e-4],np.arange(120e-4,20*ecu,20e-4))
        
        npts = len(df)
        
        if P > 0:
            for k in range(npts):
                f1 = interp1d(ecun,fcun); temp1 = np.sum(f1(df[0]*np.ones(len(yl)))*conclay[:,1])
                f2 = interp1d(ec,fc); temp2 = np.sum(f2(df[0]*np.ones(len(yl)))*conclay[:,2])
                f3 = interp1d(es,fs); temp3 = np.sum(Asb*f3(df[0]*np.ones(len(distld))))
                compch = temp1 + temp2 + temp3
                if compch < P:
                    df = df[1:]
        
        npts = len(df)
        
        # ===============ITERATIVE PROCESS TO FIND THE MOMENT - CURVATURE RELATION: ==============================
        
        msg = 0                                 # stop conditions
        
        curv   = [0]                            # curvatures
        mom    = [0]                            # moments
        ejen   = [0]                            # neutral axis
        DF     = [0]                            # force eqilibrium
        vniter = [0]                            # iterations
        coverstrain = [0]                          
        corestrain  = [0]                          
        steelstrain = [0]
        
        tol = tolerance*0.25*np.pi*(D**2)*fpc        # tolerance allowed
        x = [D/2]                                    # location of N.A.
        for k in range(npts):
            lostmomcontrol = max(mom)
            if mom[k]<(0.8*lostmomcontrol):
                msg = 4
                break
            
            F = 10*tol
            niter = -1                                  
            while abs(F)>tol:
                niter = niter + 1
                eec = (df[k]/x[niter])*(conclay[:,0]-(D-x[niter]))    # vector with the strains in the concrete
                ees = (df[k]/x[niter])*(distld-(D-x[niter]))          # vector with the strains in the steel
        
                fcunconf = interp1d(ecun,fcun)(eec)        # vector with stresses in the unconfined concr.           
                fcconf = interp1d(ec,fc)(eec)              # vector with stresses in the confinded concr.    
                fsteel = interp1d(es,fs)(ees)              # vector with stresses in the steel
                FUNCON = fcunconf*conclay[:,1]
                FCONF  = fcconf*conclay[:,2]
                FST    = Asb*fsteel;
                F      = np.sum(FUNCON) + np.sum(FCONF) + np.sum(FST) - P
                if F>0:
                    x.append(x[niter] - 0.05*x[niter])
                
                elif F<0:
                    x.append(x[niter] + 0.05*x[niter])
        
                if niter>itermax:
                    msg = 3
                    break
        
            cores = (df[k]/x[niter])*abs(x[niter]-dcore)
            TF = confined == unconfined
            if not TF:
                if cores >= ecu:
                    msg = 1
                    break
        
            elif TF:
                if df[k] >= ecu:
                    msg = 1
                    break
                
            if abs(ees[0]) > esu:
                msg = 2
                break
        
            ejen.append(x[niter])
            DF.append(x[niter])
            vniter.append(niter)
            temp = (sum(FUNCON*conclay[:,0]) + sum(FCONF*conclay[:,0]) + sum(FST*distld) - P*(D/2))/(10**6)
            mom.append(temp)
        
            if mom[k+1] < 0:
                mom[k+1] = -0.01*mom[k+1]
        
            curv.append(1000*df[k]/x[niter])
            coverstrain.append(df[k])
            corestrain.append(cores)
            steelstrain.append(ees[0])
            x[0] = x[niter]
            del x[1:]
            if msg!=0:
                break
            
        Agross = 0.25*np.pi*(D**2)
        AsLong = nbl*Asb
        LongSteelRatio  = nbl*Asb/Agross
        TransvSteelRatio = np.pi*Dh*Dh/(s*Dsp)
        AxialRatio  = P/(fpc*Agross)
        
        Mn     = interp1d(coverstrain,mom)(0.004)
        esaux    = interp1d(mom,steelstrain)(Mn)
        if esaux <- 0.015:
            Mn    = interp1d(steelstrain,mom)(-0.015)
            
        cMn = interp1d(mom,ejen)(Mn)
        
        fycurv = interp1d(steelstrain,curv)(-fy/Es)          # curvature for first yield
        fyM    = interp1d(curv,mom)(fycurv)                  # moment for first yield
        
        eqcurv = max((Mn/fyM)*fycurv,fycurv)
        
        curvbilin = np.append(np.append(0,eqcurv), curv[-1])
        mombilin  = np.append(np.append(0,Mn), mom[-1])
        
        SectionCurvatureDuctility = curv[-1]/eqcurv
        
        self.figures['3'] = {'curvbilin': curvbilin, 'mombilin': mombilin, 'curv':curv, 'mom': mom}
        
        # ===============FIND THE FORCE - DISPLACEMENT RELATION: ==============================
        
        Lsp = np.zeros(len(steelstrain))
        for j in range(len(steelstrain)):
            ffss = (-steelstrain[j]*Es)
            if ffss > fy:
                ffss = fy
            Lsp[j] = kLsp*ffss*Dbl     # Strain penetration length      
            
        kkk = min(0.2*(fsu/fy-1),0.08)
        
        if bending == 'single':
            Lp = max(kkk*L + kLsp*fy*Dbl,2*kLsp*fy*Dbl)           # Plastic hinge length
            LBE = L
        elif bending == 'double':
            Lp = max(kkk*L + kLsp*fy*Dbl,2*kLsp*fy*Dbl)           # Plastic hinge length
            LBE = L/2
        else:
            print('bending should be specified as single or double')
            sys.exit()
        
        
        # Lp = 467.5;
        # Lsp(:) = 0.5*Lp;
        
        # Moyer - Kowalsky Buckling model
        bucritMK = 0
        CuDu   = curv/eqcurv
        
        steelstrain = np.asarray(steelstrain)
        if SectionCurvatureDuctility > 4:
            
            esgr4  = -0.5*interp1d(CuDu,steelstrain)(4)   # resdidual growth strain at ductility 4
            escc   = 3*((s/Dbl)**(-2.5))                  # allowable steel compression strain
            
            esgr = np.zeros(len(steelstrain))
            for i in range(len(steelstrain)):
                
                if CuDu[i] < 1:
                    esgr[i] = 0
                    
                elif CuDu[i] < 4 and CuDu[i] > 1:
                    esgr[i] = (esgr4/4)*CuDu[i]
        
                elif CuDu[i] > 4:
                    esgr[i] = -0.5*steelstrain[i]
        
            esfl = escc-esgr
            
            self.figures['4'] = {'opt':0,'CuDu':CuDu,'steelstrain':steelstrain, 'esfl':esfl}
            if -steelstrain[-1] >= esfl[-1]:
                bucritMK = 1
                fail     = esfl + steelstrain
                failCuDuMK = interp1d(fail,CuDu)(0)
                failesfl = interp1d(fail,esfl)(0)
                failss   = -interp1d(fail,steelstrain)(0)
                self.figures['4'] = {'opt':1,'CuDu':CuDu,'steelstrain':steelstrain, 'esfl':esfl, 'failCuDuMK':failCuDuMK, 'failss': failss}
                
        # Berry - Eberhard Buckling model
        
        bucritBE = 0
        
        if AxialRatio >= 0.30:
            C0=0.006; C1=7.190; C2=3.129; C3=0.651; C4=0.227;        # model constants
        else:
            C0=0.0010; C1=7.30; C2=1.30; C3=1.30; C4=3.00;           # model constants
        
        # effective confinement ratio
        roeff = TransvSteelRatio*fyh/fpc                             
        
        # plastic rotation at the onset of bar buckling
        rotb  = C0*(1+C1*roeff)*((1+C2*P/(Agross*fpc))**(-1))*(1+C3*LBE/D+C4*Dbl*fy/D)    
        plrot = (curv-fycurv)*(Lp)/1000
        
        self.figures['5'] = {'opt':0,'CuDu':CuDu,'plrot':plrot,'rotb':rotb}
        if max(plrot) > rotb:
            bucritBE = 1
            failBE = plrot - rotb
            failplrot  = interp1d(failBE,plrot)(0)
            failCuDuBE = interp1d(failBE,CuDu)(0)
            self.figures['5'] = {'opt':1,'CuDu':CuDu,'rotb':rotb,'plrot':plrot,'failCuDuBE':failCuDuBE,'failplrot':failplrot}

        # Flexure deflection:
        displf = np.zeros(len(curv))
        mom = np.asarray(mom)
        if bending == 'single':
            for i in range(len(curv)):
                if coverstrain[i] < eccr:
                    displf[i] = curv[i]*((L/1000)**2)/3
                    
                if coverstrain[i] > eccr and curv[i] < fycurv:
                    displf[i] = curv[i] * (((L+Lsp[i])/1000)**2)/3
        
                if curv[i] >= fycurv:
                    displf[i] = (curv[i]-fycurv*(mom[i]/fyM))*(Lp/1000)*((L+Lsp[i]-0.5*Lp)/1000) + \
                                (fycurv*(((L+Lsp[i])/1000)**2)/3)*(mom[i]/fyM)
        
            Force = mom/(L/1000)
            
        elif bending == 'double':
            for i in range(len(curv)):
                if coverstrain[i] < eccr:
                    displf[i] = curv[i]*((L/1000)**2)/6
        
                if coverstrain[i] > eccr and curv[i] < fycurv:
                    displf[i] = curv[i] * (((L+2*Lsp[i])/1000)**2)/6
        
                if curv[i] >= fycurv:
                    displf[i] = (curv[i]-fycurv*(mom[i]/fyM))*(Lp/1000)*((L+2*(Lsp[i]-0.5*Lp))/1000) + \
                                (fycurv*(((L+2*Lsp[i])/1000)^2)/6)*(mom[i]/fyM);
        
            Force = 2*mom/(L/1000)
        
        else:
            print('bending should be specified as single or double')
            sys.exit()
            
        # Shear deflection:
        
        G     = 0.43*Ec
        As    = 0.9*Agross
        Ig    = np.pi*(D**4)/64
        Ieff  = (Mn*1000/(Ec*(10**6)*eqcurv))*(10**12)
        
        beta  = min(0.5+20*LongSteelRatio,1)
        
        if bending == 'single':
            alpha = min(max(1,3-L/D),1.5)
            
        elif bending == 'double':
            alpha = min(max(1,3-L/(2*D)),1.5)
        
        Vc1   = 0.29*alpha*beta*0.8*(fpc**(1/2))*Agross/1000
        
        kscr  = ((0.39*TransvSteelRatio)*0.25*Es*((0.8*D/1000)**2)/(0.25+10*(0.39*TransvSteelRatio)))*1000
        
        if bending == 'single':
            ksg   = (G*As/L)/1000
            kscr  = (kscr/L)
            forcebilin = mombilin/(L/1000)
            
        elif bending == 'double':    
            ksg   = (G*As/(L/2))/1000
            kscr  = (kscr/(L/2))
            forcebilin = 2*mombilin/(L/1000)
            
        kseff = ksg*(Ieff/Ig)
        aux = (Vc1/kseff)/1000
        aux2 = 0
        momaux = mom*1
        displsh = np.zeros(len(curv))
        for i in range(len(curv)):
            if momaux[i] <= Mn and Force[i] < Vc1:
                displsh[i] = (Force[i]/kseff)/1000
        
            if momaux[i] <= Mn and Force[i] >= Vc1:
                displsh[i] = ((Force[i]-Vc1)/kscr)/1000+aux
        
            if momaux[i] > Mn:
                momaux = 4*momaux
                aux3 = i - aux2
                aux2 = aux2 + 1
                displsh[i] = (displf[i]/displf[i-1])*displsh[i-1]
        
        displ = displsh + displf
        
        # bilinear approx:
        dy1 = interp1d(curv,displ)(fycurv)
        dy  = (Mn/fyM)*dy1
        du  = displ[-1]
        displbilin  = np.append(np.append(0,dy),du)
        Dduct = displ/dy
        DisplDuct = max(Dduct)
        dy1f = interp1d(curv,displf)(fycurv)
        dyf  = (Mn/fyM)*dy1f
        
        # Shear Strength:
        Vs    = (0.5*np.pi*(0.25*np.pi*(Dh**2))*fyh/np.tan(np.pi/6)*(D-clb+0.5*Dh-cMn)/s)/1000
        Vsd   = (0.5*np.pi*(0.25*np.pi*(Dh**2))*fyh/np.tan((35/180)*np.pi)*(D-clb+0.5*Dh-cMn)/s)/1000
        beta  = min(0.5+20*LongSteelRatio,1)
        Dductf = displ/dyf
        
        if bending == 'single':
            alpha = min(max(1,3-L/D),1.5)
            if P > 0:
                Vp = (P*(D-cMn)/(2*L))/1000
            else:
                Vp = 0
        
        elif bending == 'double':
            alpha = min(max(1,3-L/(2*D)),1.5)
            if P > 0:
                Vp = (P*(D-cMn)/(L))/100
            else:
                Vp = 0
        
        Vc = np.zeros(len(Dductf))
        if ductilitymode == 'uniaxial':
            for i in range(len(Dductf)):
                Vc[i] = alpha*beta*min(max(0.05,0.37-0.04*Dductf[i]),0.29)*0.8*(fpc**(1/2))*Agross/1000
        
        elif ductilitymode == 'biaxial':        
            for i in range(len(Dductf)):
                Vc[i] = alpha*beta*min(max(0.05,0.33-0.04*Dductf[i]),0.29)*0.8*(fpc**(1/2))*Agross/1000
                
        Vcd = 0.862*Vc
        Vpd = 0.85*Vp
        V   = Vc + Vs + Vp
        Vd  = 0.85 * (Vcd + Vsd + Vpd)
        criteria = 1
        if V[-1] < Force[-1]:
            failure   = V-Force
            faildispl = interp1d(failure,displ)(0)
            failforce = interp1d(displ,Force)(faildispl)
            failduct  = interp1d(displ,Dduct)(faildispl)
            failmom   = interp1d(displ,mom)(faildispl)
            failcurv  = interp1d(displ,curv)(faildispl)
            failCuDu  = interp1d(displ,CuDu)(faildispl)
            
            if bending == 'single':
                if faildispl <= 2*dy:
                    criteria = 2
                elif faildispl < 8*dy:
                    criteria = 3
                else:
                    criteria = 4
            elif bending == 'double':
                if faildispl <= 1*dy:
                    criteria = 2
                elif faildispl < 7*dy:
                    criteria = 3
                else:
                    criteria = 4
                    
        Ieq = (Mn/(eqcurv*Ec))/1000 # equivalent I for NLTHA
        Bi = 1/(((mombilin[1])/(curvbilin[1]))/((mombilin[2]-mombilin[1])/(curvbilin[2]-curvbilin[1]))) # Bilinear factor
        
        # Limit States:
        displdam = 0; displser = 0; Dductdam = 0; Dductser = 0;
        curvdam = 0;   curvser  = 0; CuDudam  = 0; CuDuser = 0;
        coverstraindam = 0; coverstrainser = 0;
        steelstraindam = 0; steelstrainser = 0;
        momdam = 0; momser = 0; Forcedam = 0; Forceser = 0; 
        
        if max(coverstrain) > ecser or max(abs(steelstrain)) > abs(esser):
            
            if max(coverstrain) > ecdam or max(abs(steelstrain)) > abs(esdam):
        
                displdamc = interp1d(coverstrain,displ,fill_value="extrapolate")(ecdam)
                displdams = interp1d(steelstrain,displ,fill_value="extrapolate")(esdam)
                displdam  = min (displdamc,displdams)
                Dductdam  = interp1d(displ,Dduct,fill_value="extrapolate")(displdam)
                curvdam   = interp1d(displ,curv,fill_value="extrapolate")(displdam)
                CuDudam   = interp1d(displ,CuDu,fill_value="extrapolate")(displdam)
                coverstraindam = interp1d(displ,coverstrain,fill_value="extrapolate")(displdam)
                steelstraindam = interp1d(displ,steelstrain,fill_value="extrapolate")(displdam)
                momdam   = interp1d(displ,mom,fill_value="extrapolate")(displdam)
                Forcedam = interp1d(displ,Force,fill_value="extrapolate")(displdam)
        
            displserc = interp1d(coverstrain,displ,fill_value="extrapolate")(ecser)
            displsers = interp1d(steelstrain,displ,fill_value="extrapolate")(esser)
            displser  = min(displserc,displsers)
            Dductser  = interp1d(displ,Dduct,fill_value="extrapolate")(displser)
            curvser   = interp1d(displ,curv,fill_value="extrapolate")(displser)
            CuDuser   = interp1d(displ,CuDu,fill_value="extrapolate")(displser)
            coverstrainser = interp1d(displ,coverstrain,fill_value="extrapolate")(displser)
            steelstrainser = interp1d(displ,steelstrain,fill_value="extrapolate")(displser)
            momser   = interp1d(displ,mom,fill_value="extrapolate")(displser)
            Forceser = interp1d(displ,Force,fill_value="extrapolate")(displser)
            
        outputlimit = [coverstrainser,steelstrainser,momser,Forceser,curvser,CuDuser,displser,Dductser, \
                        coverstraindam,steelstraindam, momdam, Forcedam, curvdam, CuDudam, displdam, Dductdam, \
                        max(coverstrain), min(steelstrain), mom[-1], Force[-1], max(curv), max(CuDu), max(displ), max(Dduct)]
        for i in range(len(outputlimit)):
            outputlimit[i] = float(outputlimit[i])
        
        # Plot the capacity curves
        pointsdam = np.where(displ<=displdam)[0]
        pointsser = np.where(displ<=displser)[0]
        self.figures['6'] = {'opt':0,'displ':displ,'Force':Force,'displbilin':displbilin, 'forcebilin':forcebilin, \
                             'V':V, 'Vd':Vd, 'pointsdam':pointsdam, 'pointsser':pointsser}
            
        if criteria != 1 and (bucritMK == 1 and bucritBE==0): 
            buckldispl = interp1d(CuDu,displ)(failCuDuMK)
            bucklforce = interp1d(CuDu,Force)(failCuDuMK)
    
            self.figures['6'] = {'opt':1,'displ':displ,'Force':Force,'displbilin':displbilin, 'forcebilin':forcebilin, \
                         'V':V, 'Vd':Vd, 'pointsdam':pointsdam, 'pointsser':pointsser, \
                             'faildispl': faildispl, 'failforce':failforce, 'buckldispl':buckldispl, 'bucklforce':bucklforce}

        elif criteria != 1 and (bucritMK == 1 and bucritBE==1): 
            buckldispl = interp1d(CuDu,displ)(failCuDuMK)
            bucklforce = interp1d(CuDu,Force)(failCuDuMK)
            buckldisplBE = interp1d(CuDu,displ)(failCuDuBE)
            bucklforceBE = interp1d(CuDu,Force)(failCuDuBE)

            self.figures['6'] = {'opt':2,'displ':displ,'Force':Force,'displbilin':displbilin, 'forcebilin':forcebilin, \
                         'V':V, 'Vd':Vd, 'pointsdam':pointsdam, 'pointsser':pointsser, \
                             'faildispl': faildispl, 'failforce':failforce, 'buckldispl':buckldispl, 'bucklforce':bucklforce, \
                                 'buckldisplBE':buckldisplBE, 'bucklforceBE':bucklforceBE}
            
        elif criteria !=1 and (bucritMK == 0 and bucritBE==1): 
            buckldisplBE = interp1d(CuDu,displ)(failCuDuBE)

            self.figures['6'] = {'opt':3,'displ':displ,'Force':Force,'displbilin':displbilin, 'forcebilin':forcebilin, \
                         'V':V, 'Vd':Vd, 'pointsdam':pointsdam, 'pointsser':pointsser, \
                             'faildispl': faildispl, 'failforce':failforce,'buckldisplBE':buckldisplBE, 'bucklforceBE':bucklforceBE}
        
        elif criteria !=1 and (bucritMK == 0 and bucritBE==0):
            self.figures['6'] = {'opt':4,'displ':displ,'Force':Force,'displbilin':displbilin, 'forcebilin':forcebilin, \
                         'V':V, 'Vd':Vd, 'pointsdam':pointsdam, 'pointsser':pointsser, \
                             'faildispl': faildispl, 'failforce':failforce}
        
        elif criteria == 1 and (bucritMK == 1 and bucritBE==0):
            buckldispl = interp1d(CuDu,displ)(failCuDuMK)
            bucklforce = interp1d(CuDu,Force)(failCuDuMK)

            self.figures['6'] = {'opt':5,'displ':displ,'Force':Force,'displbilin':displbilin, 'forcebilin':forcebilin, \
                         'V':V, 'Vd':Vd, 'pointsdam':pointsdam, 'pointsser':pointsser, \
                         'buckldispl':buckldispl, 'bucklforce':bucklforce}
        
        elif criteria == 1 and (bucritMK == 1 and bucritBE==1): 
            buckldispl = interp1d(CuDu,displ)(failCuDuMK)
            bucklforce = interp1d(CuDu,Force)(failCuDuMK)
            buckldisplBE = interp1d(CuDu,displ)(failCuDuBE)
            bucklforceBE = interp1d(CuDu,Force)(failCuDuBE)

            self.figures['6'] = {'opt':6,'displ':displ,'Force':Force,'displbilin':displbilin, 'forcebilin':forcebilin, \
                         'V':V, 'Vd':Vd, 'pointsdam':pointsdam, 'pointsser':pointsser, \
                         'buckldispl':buckldispl, 'bucklforce':bucklforce, 'buckldisplBE':buckldisplBE, 'bucklforceBE':bucklforceBE}
        
        elif criteria == 1 and (bucritMK == 0 and bucritBE==1): 
            buckldisplBE = interp1d(CuDu,displ)(failCuDuBE)
            bucklforceBE = interp1d(CuDu,Force)(failCuDuBE)

            self.figures['6'] = {'opt':7,'displ':displ,'Force':Force,'displbilin':displbilin, 'forcebilin':forcebilin, \
                         'V':V, 'Vd':Vd, 'pointsdam':pointsdam, 'pointsser':pointsser, \
                         'buckldisplBE':buckldisplBE, 'bucklforceBE':bucklforceBE}
        
        elif criteria == 1 and (bucritMK == 0 and bucritBE==0): 
            pass
        
        # ==========================================================================
        
        output = np.asarray([np.asarray(coverstrain), np.asarray(corestrain), np.asarray(ejen), steelstrain, \
                  mom, np.asarray(curv), Force, displsh, displf, displ, V, Vd]).T
            
        outputbilin = np.asarray([curvbilin, mombilin, displbilin, forcebilin]).T
        
        Acore = 0.25*np.pi*(Dsp**2)
        
        # compression force for yield surface
        PCid = interp1d(ec,fc)(csid)*(Acore-AsLong)+interp1d(ecun,fcun)(csid)*(Agross-Acore)+AsLong*interp1d(es,fs)(csid)
        # tensile force for yield surface
        PTid = AsLong*interp1d(es,fs)(ssid)
        
        fid = open(name +'.xlsx','w')
        
        fid.write('Circular Section\n\n')
        if confined == 'mclw':
            fid.write('Lightweight concrete\n')
        else:
            fid.write('Normalweight concrete\n')
        
        fid.write('Diameter:\t%5.1f\tmm\n' % D)
        fid.write('Cover to longitudinal bars:\t%.1f\tmm\n' % clb)
        fid.write('Number of longitudinal bars:\t%.0f\t\n' % nbl)
        fid.write('Diameter of longitudinal bars:\t%.1f\tmm\n' % Dbl)
        fid.write('Diameter of transverse steel:\t%.1f\tmm\n' % Dh)
        fid.write('Spacing of transverse steel:\t%.1f\tmm\n' % s)
        
        if stype == 'spirals':
            fid.write('Type of tranverse reinforcement: Spirals\n')
        elif stype == 'hoops':
            fid.write('Type of tranverse reinforcement: Hoops\n')
        
        fid.write('Axial load:\t%8.2f\tkN\n' % float(P/1000))
        fid.write('Concrete compressive strength:\t%3.2f\tMPa\n' % fpc)
        fid.write('Long steel yielding stress:\t%4.2f\tMPa\n' % fy)
        fid.write('Long steel max. stress:\t%4.2f\tMPa\n' % max(fs))
        fid.write('Transverse steel yielding stress:\t%4.2f\tMPa\n' % fyh)
        fid.write('Member Length:\t%5.1f\tmm\n' % L)
        
        if bending == 'single':
            fid.write('Single Bending\n')
        elif bending == 'double':
            fid.write('Double Bending\n')
        
        if ductilitymode == 'uniaxial':
            fid.write('Uniaxial Bending\n')    
        elif ductilitymode == 'biaxial':
            fid.write('Biaxial Bending\n')        
        
        fid.write('Longitudinal Steel Ratio:\t%1.3f\n' % LongSteelRatio)
        fid.write('Transverse Steel Ratio:\t%1.3f\n' % TransvSteelRatio)
        fid.write('Axial Load Ratio:\t%1.3f\n\n' % AxialRatio)
        
        fid.write('Cover\tCore\tN.A\tSteel\tMoment\tCurvature\tForce\tSh displ.\tFl displ.\tTotal displ.\tShear(assess.)\tShear(design)\n')
        fid.write('Strain\tStrain\t[mm]\tStrain\t[kN.m]\t[1/m]\t[kN]\t[m]\t[m]\t[m]\t[kN]\t[kN]\n')
        for i in range(output.shape[0]):
            fid.write('%1.5f\t%1.5f\t%4.2f\t%1.5f\t%8.2f\t%1.5f\t%8.2f\t%1.5f\t%1.5f\t%1.5f\t%8.2f\t%8.2f\n' % tuple(output[i,:]))
        fid.write('\n')
        fid.write('Bilinear Approximation:\n')
        fid.write('Curvature\tMoment\tDispl.\tForce\n')
        fid.write('[1/m]\t[kN.m]\t[m]\t[kN]\n')
        for i in range(outputbilin.shape[0]):
            fid.write('%1.5f\t%8.2f\t%1.5f\t%8.2f\n' % tuple(outputbilin[i,:]))
        fid.write('\n')
        
        if msg == 1:
           fid.write('*** Concrete strain exceeds maximum ***')
        elif msg == 2:
           fid.write('*** Steel strain exceeds maximum ***')
        elif msg == 3:
           fid.write('*** Number of iteration exceeds maximum ***')
        elif msg == 4:
           fid.write('*** Excessive loss of strength ***')
        
        fid.write('\nMoment for First Yielding:\t%8.2f\tkN.m\n' % fyM)
        fid.write('Curvature for First Yielding:\t%1.5f\t1/m\n' % fycurv)
        fid.write('Potential Section Nominal Moment:\t%8.2f\tkN.m\n' % Mn)
        fid.write('Equivalent Curvature:\t%1.5f\t1/m\n' % eqcurv)
        fid.write('Potential Section Curvature Ductility:\t%3.2f\n' % SectionCurvatureDuctility)
        fid.write('Potential Displacement Ductility:\t%3.2f\n' % DisplDuct)
        fid.write('\n')
        
        if criteria == 1:
            fid.write('*** flexural failure ***')
        elif criteria == 2:
            fid.write('*** brittle shear failure ***')
            fid.write('\nDisplacement for Shear Failure:\t%1.5f\tm\n' % faildispl)
            fid.write('Displacement Ductility at Shear Failure:\t%8.2f\n' % failduct)
            fid.write('Force for Shear Failure:\t%8.2f\tkN\n' % failforce)
            fid.write('Curvature for Shear Failure:\t%1.5f\t1/m\n' % failcurv)
            fid.write('Curvature Ductility at Shear Failure:\t%8.2f\n' % failCuDu)
            fid.write('Moment for Shear Failure:\t%8.2f\tkN.m\n' % failmom)
                
        elif criteria == 3:
            fid.write('*** shear failure at some ductility ***')
            fid.write('\nDisplacement for Shear Failure:\t%1.5f\tm\n' % faildispl);
            fid.write('Displacement Ductility at Shear Failure:\t%8.2f\n' % failduct)
            fid.write('Force for Shear Failure:\t%8.2f\tkN\n' % failforce)
            fid.write('Curvature for Shear Failure:\t%1.5f\t1/m\n' % failcurv)
            fid.write('Curvature Ductility at Shear Failure:\t%8.2f\n' % failCuDu)
            fid.write('Moment for Shear Failure:\t%8.2f\tkN.m\n' % failmom)
                
        elif criteria == 4:
            fid.write('*** ductil shear failure ***')
            fid.write('\nDisplacement for Shear Failure:\t%1.5f\tm\n' % faildispl)
            fid.write('Displacement Ductility at Shear Failure:\t%8.2f\n' % failduct)
            fid.write('Force for Shear Failure:\t%8.2f\tkN\n' % failforce)
            fid.write('Curvature for Shear Failure:\t%1.5f\t1/m\n' % failcurv)
            fid.write('Curvature Ductility at Shear Failure:\t%8.2f\n' % failCuDu)
            fid.write('Moment for Shear Failure:\t%8.2f\tkN.m\n' % failmom)
                
        if bucritMK == 1:
            bucklDd =   interp1d(CuDu,Dduct)(failCuDuMK)
            bucklcurv = interp1d(CuDu,curv)(failCuDuMK)
            bucklmom  = interp1d(CuDu,mom)(failCuDuMK)
            fid.write('Moyer - Kowalsky buckling model:\n')
            fid.write('\nCurvature Ductility for Buckling:\t%8.2f\n' % failCuDuMK)
            fid.write('Curvature at Buckling:\t%3.5f\tm\n' % bucklcurv)
            fid.write('Displacement Ductility at Buckling:\t%8.2f\n' % bucklDd)
            fid.write('Displacement at Buckling:\t%3.5f\tm\n' % buckldispl)
            fid.write('Force for Buckling:\t%8.2f\tkN\n' % bucklforce)
            fid.write('Moment for Buckling:\t%8.2f\tkN\n' % bucklmom)
            
        if bucritBE == 1:
            bucklDdBE = interp1d(CuDu,Dduct)(failCuDuBE)
            bucklcurvBE = interp1d(CuDu,curv)(failCuDuBE)
            bucklmomBE  = interp1d(CuDu,mom)(failCuDuBE)
            fid.write('Berry - Eberhard buckling model:\n')
            fid.write('\nCurvature Ductility for Buckling:\t%8.2f\n' % failCuDuBE)
            fid.write('Curvature at Buckling:\t%3.5f\tm\n' % bucklcurvBE) 
            fid.write('Displacement Ductility at Buckling:\t%8.2f\n' % bucklDdBE)
            fid.write('Displacement at Buckling:\t%3.5f\tm\n' % buckldisplBE)
            fid.write('Force for Buckling:\t%8.2f\tkN\n' % bucklforceBE)
            fid.write('Moment for Buckling:\t%8.2f\tkN\n' % bucklmomBE)
            
        fid.write('\n')
        fid.write('*** Potential Deformation Limit States (serviceability/damage control/ultimate) ***\n')
        fid.write('Cover\tSteel\tMoment\tForce\tCurvature\tCurvature\tDisplacement\tDisplacement\n')
        fid.write('Strain\tStrain\t[kN.m]\t[kN]\t[1/m]\tDuctility\t[m]\tDuctility\n')
        for i in range(1,4):
            fid.write('%1.5f\t%1.5f\t%8.2f\t%8.2f\t%1.5f\t%3.2f\t%2.5f\t%3.2f\n' % tuple(outputlimit[(i-1)*8:(8*i)]))
        fid.write('\nDeformation Limit States Citeria :\n')
        fid.write('serviceability concrete strain:\t%1.4f\n' % ecser)
        fid.write('serviceability steel strain:\t%1.4f\n' % esser)
        fid.write('damage control concrete strain:\t%1.4f\n' % ecdam)
        fid.write('damage control steel strain:\t%1.4f\n' % esdam)
        
        if confined == 'mc':
            fid.write('\nOriginal Mander Model Ultimate Concrete Strain:\t%1.4f\n' % ecumander)
        
        fid.write('\nfor NLTHA:\n')
        fid.write('E: \t%10.2f \tPa\n' % float(Ec*(10**6)))
        fid.write('G: \t%10.2f \tPa\n' % float(G*(10**6)))
        fid.write('A: \t%10.4f \tm2\n' % float(Agross/(10**6)))
        fid.write('I: \t%10.6f \tm4\n' % float(Ieq))
        fid.write('Bi-Factor: \t%1.3f\n' % float(Bi))
        fid.write('Hinge Length: \t%1.3f \tm\n' % float(Lp/1000))
        fid.write('Tension Yield: \t%10.2f \tN\n' % float(PTid))
        fid.write('Compression Yield: \t%10.2f \tN\n' % float(PCid))
        fid.write('Moment Yield: \t%10.2f \tN-m\n' % float(Mn*1000))
        
        # vector with axial loads for interaction diagram
        temp1 = np.append(np.arange(-0.90*PTid,0,0.30*PTid),0).tolist()
        temp2 = np.append(np.arange(0.05*fpc*Agross,0.6*PCid,0.05*fpc*Agross),0.6*PCid).tolist()
        PP = temp1 + temp2 +[0.7*PCid] + [0.8*PCid] + [0.9*PCid]
        nPP = len(PP)
        Mni = []
        eci = []
        esi = []
        Msgs = []
        
        for i in range(nPP):
            
            df = np.arange(1e-4,20*ecu,1e-4)
                
            if ecu > 0.0018:
                df = np.append(df[df<=16e-4],np.arange(18e-4,20*ecu,2e-4))
                
            if ecu > 0.0025:
                df = np.append(df[df<=20e-4],np.arange(25e-4,20*ecu,5e-4))
                
            if ecu > 0.006:
                df = np.append(df[df<=50e-4],np.arange(60e-4,20*ecu,10e-4))
                
            if ecu > 0.012: 
                df = np.append(df[df<=100e-4],np.arange(120e-4,20*ecu,20e-4))
            
            npts = len(df)
            
            if PP[i] > 0:
                for k in range(npts):
                    f1 = interp1d(ecun,fcun); temp1 = np.sum(f1(df[0]*np.ones(len(yl)))*conclay[:,1])
                    f2 = interp1d(ec,fc); temp2 = np.sum(f2(df[0]*np.ones(len(yl)))*conclay[:,2])
                    f3 = interp1d(es,fs); temp3 = np.sum(Asb*f3(df[0]*np.ones(len(distld))))
                    compch = temp1 + temp2 + temp3
                    if compch < PP[i]:
                        df = df[1:]
            
            npts = len(df)
        
            msg = 0                                 # stop conditions
            
            curv   = [0]                            # curvatures
            mom    = [0]                            # moments
            ejen   = [0]                            # neutral axis
            DF     = [0]                            # force eqilibrium
            vniter = [0]                            # iterations
            coverstrain = [0]                          
            corestrain  = [0]                          
            steelstrain = [0]
            
            x = [D/2]
            for k in range(npts):
                lostmomcontrol = max(mom)
                if mom[k]<(0.8*lostmomcontrol):
                    msg = 4
                    break
                
                F = 10*tol
                niter = -1                                  
                while abs(F)>tol:
                    niter = niter + 1
                    eec = (df[k]/x[niter])*(conclay[:,0]-(D-x[niter]))    # vector with the strains in the concrete
                    ees = (df[k]/x[niter])*(distld-(D-x[niter]))          # vector with the strains in the steel
            
                    fcunconf = interp1d(ecun,fcun,fill_value="extrapolate")(eec)        # vector with stresses in the unconfined concr.           
                    fcconf = interp1d(ec,fc,fill_value="extrapolate")(eec)              # vector with stresses in the confinded concr.    
                    fsteel = interp1d(es,fs,fill_value="extrapolate")(ees)              # vector with stresses in the steel
                    FUNCON = fcunconf*conclay[:,1]
                    FCONF  = fcconf*conclay[:,2]
                    FST    = Asb*fsteel
                    F      = np.sum(FUNCON) + np.sum(FCONF) + np.sum(FST) - PP[i]
                    
                    if F>0:
                        x.append(x[niter] - 0.05*x[niter])
                    
                    elif F<0:
                        x.append(x[niter] + 0.05*x[niter])
            
                    if niter>itermax:
                        msg = 3
                        break
            
                cores = (df[k]/x[niter])*abs(x[niter]-dcore)
                TF = confined == unconfined
                if not TF:
                    if cores >= ecu:
                        msg = 1
                        break
            
                elif TF:
                    if df[k] >= ecu:
                        msg = 1
                        break
                    
                if abs(ees[0]) > esu:
                    message = 2
                    break
            
                ejen.append(x[niter])
                DF.append(x[niter])
                vniter.append(niter)
                temp = (sum(FUNCON*conclay[:,0]) + sum(FCONF*conclay[:,0]) + sum(FST*distld) - PP[i]*(D/2))/(10**6)
                mom.append(temp)
            
                if mom[k+1] < 0:
                    mom[k+1] = -0.01*mom[k+1]
            
                curv.append(1000*df[k]/x[niter])
                coverstrain.append(df[k])
                corestrain.append(cores)
                steelstrain.append(ees[0])
                x[0] = x[niter]
                del x[1:]
                if msg!=0:
                    break
            
            Mni.append(interp1d(coverstrain,mom)(csid))
            esaux    = interp1d(coverstrain,steelstrain)(csid)
            cr = 0 # concrete control
            if abs(esaux) > abs(ssid) or np.isnan(Mni[i]):
                cr = 1 # steel concrete
                Mni[i]    = interp1d(steelstrain,mom)(-ssid)
            
            if cr == 0:
                eci.append(csid)
                esi.append(esaux)
            if cr == 1:
                esi.append(-ssid)
                eci.append(float(interp1d(steelstrain,coverstrain)(-ssid)))    
                
            Msgs.append(msg)
        
        Mni = np.asarray([0]+Mni+[0])
        PPn  = np.asarray([-PTid] + PP + [PCid])
        
        MB = np.max(Mni)
        PB = PPn[np.where(Mni==MB)[0]][0]
        
        PB13 = (1/3)*PB
        MB13 = float(interp1d(PPn,Mni)(PB13))
        
        PB23 = (2/3)*PB
        MB23 = float(interp1d(PPn,Mni)(PB23))
        
        MB0 = float(interp1d(PPn,Mni)(0))
        
        PPL = np.asarray([-PTid, 0, PB13, PB23, PB, PCid])
        MnL = np.asarray([0, MB0, MB13, MB23, MB, 0])

        self.figures['7'] = {'Mni':Mni,'PPn':PPn,'MnL':MnL,'PPL':PPL}

        fid.write('\n\n*** Interaction Surface ***\n')
        fid.write('Concrete limit strain:\t%1.4f\n' % csid)
        fid.write('Steel limit strain:\t%1.4f\n' % ssid)
        fid.write('\nMoment\tAxial Load\n') 
        fid.write('[kN-m]\t[kN]\n')
        for i in range(len(Mni)):
            fid.write('%8.2f\t%8.2f\n' % (Mni[i], PPn[i]))
            
        fid.write('\n')
        fid.write('NLTHA Approximation:\n\n')
        fid.write('PT:\t%6.1f\tkN\n' % float(-PTid/1000))
        fid.write('PC:\t%6.1f\tkN\n' % float(PCid/1000))
        fid.write('PB:\t%6.1f\tkN\tMB:\t%6.1f\tkN-m\n' % (float(PB/1000),MB))
        fid.write('(1/3)PB:\t%6.1f\tkN\t(1/3)MB:\t%6.1f\tkN.m\n' % (float(PB13/1000),MB13))
        fid.write('(2/3)PB:\t%6.1f\tkN\t(2/3)MB:\t%6.1f\tkN.m\n' % (float(PB23/1000),MB23))

        self.figures['8'] = {'PPn':PPn,'eci':eci}
        
        self.figures['9'] = {'PPn':PPn,'esi':eci}
        
        fid.close()

''' Activate '''
app = Gui()
app.root.mainloop()