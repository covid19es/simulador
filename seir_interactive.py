import math
import numpy as np
import pandas as pd
from bokeh.plotting import figure, output_file, show
from bokeh.models import ColumnDataSource
from bokeh.models.tools import HoverTool
from bokeh.io import output_notebook
from bokeh.models import Legend
from bokeh.palettes import Spectral4
from bokeh.layouts import column, row, WidgetBox
from bokeh.models import Panel
from bokeh.models.widgets import Tabs
from bokeh.models.widgets import Slider


from bokeh.application.handlers import FunctionHandler
from bokeh.application import Application

output_notebook() 

Time_default              = 200 # número de días de simulación

InterventionTime_default  = 100 #dias hasta el confinamiento 
InterventionAmt_default   = 4/5 # efectividad del confinamiento (valor de la reducción)


Ro_default                = 5.55 # tasa inicial de infección del virus

D_incbation_default      = 5.2  #tiempo de incubación      
D_infectious_default      = 3.2 #tiempo que un paciente infecta

CFR_array_default = np.array([0.0016, 0.0045, 0.0892]) # valor original
CFR_array_default = np.array([0.002, 0.005, 0.089]) # valor redondeado

P_SEVERE_array_default = np.array([0.0098, 0.0179, 0.0528]) # valor original
P_SEVERE_array_default = np.array([0.001, 0.017, 0.053]) # valor redondeado


#criterios clínicos
Time_to_death_default     = 32 - D_incbation_default #numero de días desde que acaba el periodo de incubación hasta la muerta
D_hospital_lag_default    = 5 # tiempo hasta la hospitalización
# NB: ambos siguientes sin ajustar por tiempo de incubación, movido a make_dataset
D_recovery_mild_default   = 14 - D_incbation_default  #tiempo en recuperarse una persona con síntomas normales
D_recovery_severe_default = 31.5 - D_incbation_default #tiempo en recuperarse una persona con síntomas severos


I0                = 1 # número de infectados
dt                = 1 #número de días por timestep
n_groups          = 3 # número de grupos poblacionales a considerar

N_young   = [2038663, 2355494, 2491371, 2322446]  # 0-4 años, 5-9 años, 10-14 años, 15-19 años
N_adult   = [2305303, 2566361, 2860390, 3492600, 4013566, 3818921, 3637840, 3266949]  # 20-59 años
N_elderly = [2799111, 2398795, 2175840, 1618039, 1353814, 975247, 415548, 103607, 16303]  # 60+ años

region_list = ["España", 
              "Andalucía",
              "Aragon",
              "Asturias",
              "Islas Baleares",
              "Canarias",
              "Cantabria",
              "Castilla y León",
              "Castilla-La Mancha",
              "Cataluña",
              "Comunidad Valenciana",
              "Extremadura",
              "Galicia",
              "Comunidad de Madrid",
              "Murcia",
              "Navarra",
              "País Vasco",
              "La Rioja",
              "Ceuta",
              "Melilla"]

# Fuente: INE
# Población para cada segmento poblacional por regiones (total + orden de CCAAs en region_list)
pop_array=[[9254359, 25802699, 12043406],
[1783728,4703083,1959763],
[249083,703769,371541],
[150052,526806,343135],
[235310,706899,256368],
[402930,1324194,493152],
[102990,309096,169868],
[388391,1234010,780469],
[411859,1119610,506989],
[1563305,4144209,1901986],
[990177,2723348,1285201],
[197720,576615,288454],
[428354,1406222,864189],
[1360320,3748132,1577018],
[339811,838411,316221],
[135096,348188,169253],
[406390,1131155,644378],
[60686,166994,86809],
[22505,47029,14900],
[25652,44929,13712]]

N_array_default = pop_array[0]  # población

RK4               = [[1]]

from bokeh.plotting import figure, output_file, show
from bokeh.models import ColumnDataSource
from bokeh.models.tools import HoverTool
from bokeh.io import output_notebook
from bokeh.models import Legend
from bokeh.palettes import Spectral4
from bokeh.layouts import column, row, WidgetBox
from bokeh.models import Panel
from bokeh.models import Range1d

from bokeh.models.widgets import Slider, Button, CheckboxGroup,Tabs, Select

from bokeh.application.handlers import FunctionHandler
from bokeh.application import Application

from functools import reduce

def modify_doc(doc):    
    
    # método para crear el dataset que será mostrado en la gráfica
    def make_dataset(Ro=Ro_default,
                     InterventionTime = InterventionTime_default, 
                     InterventionAmt = InterventionAmt_default,
                     Time_to_death = Time_to_death_default,
                     D_recovery_mild=D_recovery_mild_default,
                     D_recovery_severe = D_recovery_severe_default,
                     D_hospital_lag = D_hospital_lag_default,
                     D_incbation = D_incbation_default,
                     D_infectious = D_infectious_default,
                     CFR_array = CFR_array_default,
                     P_SEVERE_array = P_SEVERE_array_default,
                     N_array=N_array_default,
                     Time = Time_default,
                     edades_seleccionadas = [0,1,2]):
                        
        #función que computa un step del vector de estados
        def f(stepFloat,v,Rvector,dt2, idx):
            #obtenemos los parametros del vector de inicio
            S        = v[0] # Susectable
            E        = v[1] # Exposed
            I        = v[2] # Infectious 
            Mild     = v[3] # Recovering (Mild)     
            Severe   = v[4] # Recovering (Severe at home)
            Severe_H = v[5] # Recovering (Severe in hospital)
            Fatal    = v[6] # Recovering (Fatal)
            R_Mild   = v[7] # Recovered
            R_Severe = v[8] # Recovered
            R_Fatal  = v[9] # Dead 

            beta2 = Rvector[int(stepFloat)]/(D_infectious) 
            p_severe = P_SEVERE_array[idx]
            p_fatal  = CFR_array[idx]
            p_mild   = 1 - P_SEVERE_array[idx] - CFR_array[idx]
            dS        = -beta2*I*S
            dE        =  beta2*I*S - a*E
            dI        =  a*E - gamma*I
            dMild     =  p_mild*gamma*I   - (1/D_recovery_mild_array[idx])*Mild
            dSevere   =  p_severe*gamma*I - (1/D_hospital_lag_array[idx])*Severe
            dSevere_H =  (1/D_hospital_lag_array[idx])*Severe - (1/D_recovery_severe_array[idx])*Severe_H
            dFatal    =  p_fatal*gamma*I  - (1/D_death_array[idx])*Fatal
            dR_Mild   =  (1/D_recovery_mild_array[idx])*Mild
            dR_Severe =  (1/D_recovery_severe_array[idx])*Severe_H
            dR_Fatal  =  (1/D_death_array[idx])*Fatal   

            return np.array([dS, dE, dI, dMild, dSevere, dSevere_H, dFatal, dR_Mild, dR_Severe, dR_Fatal])

        def integrate(m,f,y,t,h,Rvector, idx):
            k =np.empty((len(m),len(y)))
            for ki in range(0,len(m)):
                _y=np.array(y.copy())
                dt=((m[ki-1][0])*h)
                for l in range(0,len(_y)):
                    for j in range(1,ki):
                        _y[l]=_y[l]+h*(m[ki-1][j])*(k[ki-1][l])
                k[ki]=f(t,_y,Rvector,dt, idx); 
            r = y. copy()
            for l in range(0,len(_y)):
                for j in range(0,len(k)):
                    r[l]=r[l]+h*(k[j][l])*(m[ki-1][j]);
            return r;
        
        # Cambiado para evitar errores en valores extremos 
        D_death = Time_to_death #días desde que acaba el periodo de infección hasta la muerte
       
        Time_to_death_array = np.full(n_groups, Time_to_death)
        D_recovery_mild_array = np.full(n_groups, D_recovery_mild)
        D_recovery_severe_array = np.full(n_groups, D_recovery_severe)
        D_hospital_lag_array = np.full(n_groups, D_hospital_lag)
        D_death_array = np.full(n_groups, D_death)
        
        Ro_array = np.full(n_groups, Ro)  # tasa inicial de infección del virus
        
        a = 1/D_incbation
        gamma = 1/D_infectious
                
        df_state_list = []
        dfConsolidated_list = []
        
        for idx in edades_seleccionadas:
            # print("Generating data for group", idx)
            steps = int(Time/dt)
            sample_step = 1
                        
            #creo el vector de R, con la intervención incluida
            # Miquel - parametros hardcodeados, cambiados a utilizar los globales otra vez
            Rvector = np.full((steps), Ro_array[idx])
            for i in range(InterventionTime, steps):
                Rvector[i]=(1-InterventionAmt)*Ro
                

            #creamos el vector de estado de inicio v (V0)
            #'susceptible','expuesto','infectado','recuperado leve','recuperando severo encasa','recuperando severo hosp',
            #'recuperando y morirá','recuperado estado normal','recuperado severo','muertos'
            v=np.array([(1 - I0/N_array[idx]), 0, I0/N_array[idx], 0, 0, 0, 0, 0, 0, 0])
            t = 0

            #inicio vacíos los DataFrame que utilizaremos
            df_state = pd.DataFrame()
            dfConsolidated = pd.DataFrame()

            init_steps = steps
            step=0
            while(steps):
                steps= steps-1

                dfConsolidatedDelta=pd.DataFrame([ np.round(np.array([N_array[idx]*v[9], N_array[idx]*(v[5]+v[6]),
                                                                      N_array[idx]*(v[7] + v[8]), N_array[idx]*v[2],
                                                                      N_array[idx]*v[1] ,int(step*dt)]))])
                dfConsolidatedDelta.columns=['DEAD','HOSPITAL','RECOVERED','INFECTED','EXPOSED','TIMESTEP']
                dfConsolidated=dfConsolidated.append(dfConsolidatedDelta)
                vDF=pd.DataFrame([np.round(np.array(v)*N_array[idx])])
                vDF.columns=['SUSCEPTIBLES','EXPUESTOS','INFECTADOS','RINGM','RINGS','RINGSH','RM','RS','RSH','MUERTOS']
                df_state=df_state.append(vDF)

                v =integrate(RK4,f,v,step,dt,Rvector, idx)
                t=t+dt

                step=step+1
            
            df_state_list.append(df_state)
            dfConsolidated_list.append(dfConsolidated)
        
        # cambio Miquel - sumar los grupos activos:
        suma = pd.DataFrame()
                
        suma = reduce(lambda x, y: x.add(y, fill_value=0), dfConsolidated_list)
        suma['TIMESTEP'] /= len(dfConsolidated_list)
        
        return ColumnDataSource(suma)

    
    # método para crear la gráfica
    def make_plot(source2):
        timesteps = source2.data['TIMESTEP'].tolist()
        timesteps=[str(i) for i in timesteps]

        tooltips = [
            ("Timestep", "@TIMESTEP"),
            ("Dead", "@DEAD"),
            ("Hospital", "@HOSPITAL"),
            ("Infected", "@INFECTED"),
            ("Exposed", "@EXPOSED"),
        ]

        p2 = figure()
        p2.vbar_stack(stackers=['DEAD','HOSPITAL','INFECTED','EXPOSED'],x='TIMESTEP',
                      source=source2,width=0.5, color=[Spectral4[3], Spectral4[0], Spectral4[2], Spectral4[1]],
                     legend_label = ['Dead','Hospital','Infected','Exposed'])
        p2.plot_height=400
        p2.plot_width=680
        #p2.y_range=Range1d(0,3e7)
        

        hover = HoverTool()
        hover.tooltips = tooltips

        p2.add_tools(hover)


        return p2
    
    #def update(attr, old, new):
    #   a = ro_slider.value
    
    def reset_params(event):
        ro_slider.value=Ro_default
        InterventionTime_slider.value=InterventionTime_default 
        InterventionAmt_slider.value=InterventionAmt_default
        Time_to_death_slider.value=Time_to_death_default
        D_recovery_mild_slider.value=D_recovery_mild_default
        D_recovery_severe_slider.value=D_recovery_severe_default
        D_hospital_lag_slider.value=D_hospital_lag_default
        D_incbation_slider.value=D_incbation_default
        D_infectious_slider.value=D_infectious_default
        CFR_array_new=CFR_array_default
        P_SEVERE_array_new=P_SEVERE_array_default
        Time_slider.value=Time_default
        edades_selection.active = [0,1,2]
        
        for i in range(3):
            CFR_slider[i].value = CFR_array_default[i]*100
            P_SEVERE_slider[i].value = P_SEVERE_array_default[i]*100
        
    # método ejecutado en el callback cuando se aprieta el boton para actualizar la gráfica
    def update_plot(event):
        index = region_list.index(ccaa_select.value)

        region_pop_array = pop_array[index]
        
                
        CFR_array_new =[]
        P_SEVERE_array_new = []
        for i in range(3):
            CFR_array_new.append(CFR_slider[i].value/100)
            P_SEVERE_array_new.append(P_SEVERE_slider[i].value/100)
            

        new_dataset = make_dataset(Ro=ro_slider.value,
                     InterventionTime = InterventionTime_slider.value, 
                     InterventionAmt = InterventionAmt_slider.value,
                     Time_to_death = Time_to_death_slider.value,
                     D_recovery_mild=D_recovery_mild_slider.value,
                     D_recovery_severe = D_recovery_severe_slider.value,
                     D_hospital_lag = D_hospital_lag_slider.value,
                     D_incbation = D_incbation_slider.value,
                     D_infectious = D_infectious_slider.value,
                     CFR_array = CFR_array_new,
                     P_SEVERE_array = P_SEVERE_array_new,
                     Time = Time_slider.value,
                     N_array=region_pop_array,
                     edades_seleccionadas = edades_selection.active)
        
        dataset.data.update(new_dataset.data)

    # crear elementos interactivos
    edades_labels=["0-19 años", "20-59 años", "59+ años"]
    
    edades_selection = CheckboxGroup(labels=edades_labels, height_policy="max",
                                  active = [0,1,2])
    
    ro_slider = Slider(start = 0.0, end = 10.0, 
                         step = 0.1, value = Ro_default,
                         title = 'R0 (número básico de reproducción)')
    
    InterventionTime_slider = Slider(start = 0.0, end = 200.0, 
                         step = 1.0, value = InterventionTime_default,
                         title = 'Dias hasta intervención')
    
    InterventionAmt_slider = Slider(start = 0.0, end = 1.0, 
                         step = 0.05, value = InterventionAmt_default,
                         title = 'Efectividad de la intervención')
        
    Time_slider = Slider(start = 50.0, end = 300.0, 
                         step = 1.00, value = Time_default,
                         title = 'Días a simular')
    
    D_incbation_slider = Slider(start = 1.2, end = 20.0, 
                         step = 0.2, value = D_incbation_default,
                         title = 'Días de incubación')
    
    D_infectious_slider = Slider(start = 1.2, end = 20.0, 
                         step = 0.2, value = D_infectious_default,
                         title = 'Días que un paciente es infeccioso')
    
    CFR_slider=[]
    P_SEVERE_slider=[] 
    
    for i in range(3):
        CFR_slider.append(Slider(start = 0.0, end =20, 
                         step = 0.1, value = CFR_array_default[i]*100,
                         title = 'Porcentaje de tasa de mortalidad (' + edades_labels[i] +')'))
        P_SEVERE_slider.append(Slider(start = 0.0, end = 20, 
                         step = 0.1, value = P_SEVERE_array_default[i]*100,
                         title = 'Porcentaje de casos severos (' + edades_labels[i] +')'))
    
    Time_to_death_slider = Slider(start = 1.2, end = 40.0, 
                         step = 0.5, value = Time_to_death_default,
                         title = 'Días desde mostrar síntomas hasta la muerte')
    
    D_recovery_mild_slider = Slider(start = 1.0, end = 40.0, 
                         step = 0.5, value = D_recovery_mild_default,
                         title = 'Días de recuperación con síntomas leves')
    
    D_recovery_severe_slider = Slider(start = 1.2, end = 40.0, 
                         step = 0.5, value = D_recovery_severe_default,
                         title = 'Días de recuperación con síntomas severos')
    
    D_hospital_lag_slider = Slider(start = 1.2, end = 40.0, 
                         step = 0.5, value = D_hospital_lag_default,
                         title = 'Tiempo hasta la hospitalización')


    controls = WidgetBox(ro_slider, 
                         D_incbation_slider,
                         D_infectious_slider,
                         Time_to_death_slider,
                         D_recovery_mild_slider,
                         D_recovery_severe_slider,
                         D_hospital_lag_slider)
    
    row3 =  row(InterventionTime_slider, 
                         InterventionAmt_slider, 
                         Time_slider)

    button_apply = Button(label="Simular", button_type="success", width_policy="min")
    button_apply.on_click(update_plot)
    
    button_reset = Button(label="Reiniciar parametros", button_type="warning", width_policy="min")
    button_reset.on_click(reset_params)
    
    ccaa_select = Select(title="Región",
                         value=region_list[0], 
                         options=region_list)
    
    row1 = row(column(edades_selection, ccaa_select, row(button_apply, button_reset)),column(CFR_slider[0], 
                                       CFR_slider[1], 
                                       CFR_slider[2],), column(P_SEVERE_slider[0], 
                                                               P_SEVERE_slider[1], 
                                                               P_SEVERE_slider[2]))
    
    # inicializar plot
    dataset = make_dataset()
    plot = make_plot(dataset)

    # dibujar controles y plot
    row2 = row(controls, plot)

    layout=column(row1,row2, row3)
    doc.add_root(layout)

handler = FunctionHandler(modify_doc)
app = Application(handler)
