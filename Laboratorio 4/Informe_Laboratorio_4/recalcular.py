"""Recalcula tablas y ajustes; usa errores estándar de la regresión."""
from pathlib import Path
import csv, json, math

ROOT = Path(__file__).resolve().parent
CONFIG = json.loads((ROOT/'datos/parametros_experimento.json').read_text(encoding='utf-8'))
with (ROOT/'datos/lecturas.csv').open(encoding='utf-8', newline='') as stream:
    RAW = list(csv.DictReader(stream))
for row in RAW:
    for key in ['f_nom','dh','time_us','dH','vdiv','probe','ymag','xmag']:
        row[key] = float(row[key])

def texnum(x, digits=6):
    return '\\num{'+format(x,f'.{digits}g')+'}'

def pair(value, uncertainty, sig=2):
    place = sig-1-math.floor(math.log10(abs(uncertainty))) if uncertainty else 3
    decimals=max(0,place)
    return (r'\ensuremath{\num{'+f'{round(value,place):.{decimals}f}'
            +r'}\pm\num{'+f'{round(uncertainty,place):.{decimals}f}'+r'}}')

def ols(x,y):
    n=len(x); xm=sum(x)/n; ym=sum(y)/n
    ss=sum((t-xm)**2 for t in x)
    a=sum((xx-xm)*(yy-ym) for xx,yy in zip(x,y))/ss
    b=ym-a*xm
    residuals=[yy-(a*xx+b) for xx,yy in zip(x,y)]
    sse=sum(r*r for r in residuals); s2=sse/(n-2)
    return dict(a=a,b=b,ua=math.sqrt(s2/ss),ub=math.sqrt(s2*(1/n+xm*xm/ss)),
                residuals=residuals,n=n)

RESULTS={}; FITS={}
for circuit in ['RC','RL']:
    R=CONFIG['R_ohm']; output=[]
    sx='VC' if circuit=='RC' else 'VL'
    for fn in sorted({r['f_nom'] for r in RAW if r['circuit']==circuit}):
        rr={r['signal']:r for r in RAW if r['circuit']==circuit and r['f_nom']==fn}
        g=rr['VG']
        T=g['dh']*g['time_us']*1e-6/g['xmag']*CONFIG['factor_tiempo']
        f=1/T
        volts={signal:z['dH']*z['vdiv']*z['probe']/z['ymag']*CONFIG['factor_vertical']/2
               for signal,z in rr.items()}
        reactance=R*volts[sx]/volts['VR']
        nominal=(1/(2*math.pi*f*CONFIG['C_nom_F']) if circuit=='RC'
                 else 2*math.pi*f*CONFIG['L_nom_H'])
        output.append(dict(circuit=circuit,f_nom=fn,dh=g['dh'],time_us=g['time_us'],
            T_us=T*1e6,f_Hz=f,f_kHz=f/1000,invf_ms=1000/f,
            VG=volts['VG'],VR=volts['VR'],**{sx:volts[sx]},X=reactance,
            X_kOhm=reactance/1000,X_nom=nominal,E_X=100*abs((reactance-nominal)/nominal)))
    RESULTS[circuit]=output
    x=[1/r['f_Hz'] if circuit=='RC' else r['f_Hz'] for r in output]
    fit=ols(x,[r['X'] for r in output])
    fit['parameter']=1/(2*math.pi*fit['a']) if circuit=='RC' else fit['a']/(2*math.pi)
    # Transformación del error estándar de la pendiente al parámetro físico.
    fit['u_fit']=fit['parameter']*fit['ua']/abs(fit['a'])
    fit['nominal']=CONFIG['C_nom_F'] if circuit=='RC' else CONFIG['L_nom_H']
    fit['E_percent']=100*abs((fit['parameter']-fit['nominal'])/fit['nominal'])
    fit['uncertainty_basis']='Error estándar por dispersión de la regresión lineal con intercepto libre.'
    FITS[circuit]=fit
    with (ROOT/f'datos/{circuit.lower()}.csv').open('w',encoding='utf-8',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(output[0]))
        writer.writeheader(); writer.writerows(output)

def table(filename,caption,label,headers,rows,note='',compact=False):
    spacing='3pt' if compact else '5pt'
    text='\\begin{center}\n\\begin{minipage}{\\linewidth}\\centering\n\\small\n\\setlength{\\tabcolsep}{'+spacing+'}\n'
    text+=r'\captionof{table}{'+caption+'}\n'+r'\label{'+label+'}\n'
    text+=r'\begin{tabular}{'+'c'*len(headers)+'}\n\\toprule\n'
    text+=' & '.join(headers)+r' \\'+'\n\\midrule\n'
    text+='\n'.join(' & '.join(row)+r' \\' for row in rows)
    text+='\n\\bottomrule\n\\end{tabular}\n'
    if note:text+='\\par\\vspace{3pt}\\begin{minipage}{.96\\linewidth}\\footnotesize '+note+'\\end{minipage}\n'
    text+='\\end{minipage}\n\\end{center}\n'
    (ROOT/'tablas'/filename).write_text(text,encoding='utf-8')

heads=[r'\shortstack{$f_{\rm nom}$\\(kHz)}',r'\shortstack{$\Delta h$\\(div)}',
       r'\shortstack{$s_t$\\($\mu$s/div)}',r'\shortstack{$T$\\($\mu$s)}',r'\shortstack{$f_{\rm exp}$\\(kHz)}']
for c in ['RC','RL']:
    freqrows=[]
    for r in RESULTS[c]:
        freqrows.append([texnum(r['f_nom']/1000),texnum(r['dh']),texnum(r['time_us']),
                         texnum(r['T_us'],4),texnum(r['f_kHz'],4)])
    table(c.lower()+'_frecuencias.tex',f'Circuito {c}: períodos y frecuencias.',f'tab:{c.lower()}_frecuencias',heads,freqrows,
          r'$s_t$: escala de tiempo por división.',compact=True)

for c in ['RC','RL']:
    symbol='C' if c=='RC' else 'L'; sx='V'+symbol
    headers=[r'\shortstack{$f_{\rm nom}$\\(kHz)}',r'\shortstack{$V_{G,m}$\\(V)}',
             r'\shortstack{$V_{R,m}$\\(V)}',r'\shortstack{$V_{'+symbol+r',m}$\\(V)}',
             r'\shortstack{$X_{'+symbol+r',\rm exp}$\\($\Omega$)}']
    rows=[[texnum(r['f_nom']/1000),texnum(r['VG']),texnum(r['VR']),texnum(r[sx]),
           texnum(r['X'],4)] for r in RESULTS[c]]
    table(c.lower()+'_resultados.tex',f'Circuito {c}: amplitudes y reactancias.','tab:'+c.lower(),headers,rows,
          r'$V_m=V_{pp}/2$.',compact=True)
    rawrows=[]
    for fn in sorted({r['f_nom'] for r in RAW if r['circuit']==c}):
        rr={r['signal']:r for r in RAW if r['circuit']==c and r['f_nom']==fn}
        g,r,v=rr['VG'],rr['VR'],rr[sx]
        rawrows.append([texnum(fn/1000),texnum(g['dH']),texnum(g['vdiv']),texnum(r['dH']),
                        texnum(r['vdiv']),texnum(v['dH']),texnum(v['vdiv'])])
    table(c.lower()+'_escalas.tex',f'Lecturas verticales de {c}. Todas las alturas son pico a pico.',
          f'tab:{c.lower()}_escalas',
          [r'$f_{\rm nom}$ (kHz)',r'$\Delta H_G$ (div)',r'$s_G$ (V/div)',
           r'$\Delta H_R$ (div)',r'$s_R$ (V/div)',r'$\Delta H_'+symbol+r'$ (div)',r'$s_'+symbol+r'$ (V/div)'],rawrows)

errors=[[texnum(rc['f_nom']/1000),texnum(rc['E_X'],3),texnum(rl['f_nom']/1000),texnum(rl['E_X'],3)]
        for rc,rl in zip(RESULTS['RC'],RESULTS['RL'])]
table('errores_reactancias.tex','Errores relativos de las reactancias en RC y RL.','tab:errores_reactancias',
      [r'\shortstack{$f_{\rm nom,RC}$\\(kHz)}',r'$E_C$ (\%)',r'\shortstack{$f_{\rm nom,RL}$\\(kHz)}',r'$E_L$ (\%)'],errors,
      r'$E_C$ corresponde a $X_C$ y $E_L$ a $X_L$.',compact=True)

fitrows=[]
for c in ['RC','RL']:
    q=FITS[c]; conv,unit=(1e6,r'$\mu$F') if c=='RC' else (1e3,'mH')
    place=-math.floor(math.log10(q['u_fit']*conv)) if q['u_fit'] else 3
    fitrows.append(['$C$' if c=='RC' else '$L$',texnum(round(q['parameter']*conv,place))+' '+unit,
                    texnum(q['u_fit']*conv,1)+' '+unit,texnum(q['E_percent'],3)])
table('ajustes.tex','Capacitancia e inductancia obtenidas de los ajustes lineales.','tab:ajustes',
      ['Magnitud','Valor obtenido',r'$u_{\rm ajuste}$',r'$E$ (\%)'],fitrows,
      r'$u_{\rm ajuste}$ es el error estándar del parámetro debido a la dispersión de los puntos en la regresión. '
      r'$E=100\left|(p_{\rm exp}-p_{\rm nom})/p_{\rm nom}\right|$.')

coeff=[[c,pair(FITS[c]['a'],FITS[c]['ua']),pair(FITS[c]['b'],FITS[c]['ub'])] for c in ['RC','RL']]
table('coeficientes.tex','Pendientes e interceptos de los ajustes de $X_C$ frente a $1/f$ y de $X_L$ frente a $f$.',
      'tab:coeficientes',
      ['Circuito',r'Pendiente $a_C$ o $a_L$',r'Intercepto $b_C$ o $b_L$ ($\Omega$)'],coeff,
      r'$a_C$ se expresa en $\Omega/\mathrm{s}$ y $a_L$ en $\Omega\,\mathrm{s}$. '
      r'Los valores después de $\pm$ son los errores estándar de los coeficientes del ajuste.')

macros={'CnomUF':CONFIG['C_nom_F']*1e6,'LnomMH':CONFIG['L_nom_H']*1e3,
        'Arcplot':FITS['RC']['a']/1e6,'Brcplot':FITS['RC']['b']/1000,
        'Arlplot':FITS['RL']['a'],'Brlplot':FITS['RL']['b']/1000}
macros_tex={k:format(v,'.12g') for k,v in macros.items()}
macros_tex.update({
    'PendienteC':format(FITS['RC']['a'],'.5g'),
    'PendienteL':format(FITS['RL']['a'],'.4g'),
    'ResultadoC':pair(FITS['RC']['parameter']*1e6,FITS['RC']['u_fit']*1e6,1),
    'ResultadoL':pair(FITS['RL']['parameter']*1e3,FITS['RL']['u_fit']*1e3,1),
})
for circuit,symbol in [('RC','C'),('RL','L')]:
    intercept=FITS[circuit]['b']
    macros_tex['Termino'+symbol]=('+' if intercept>=0 else '-')+r'\SI{'+format(abs(intercept),'.3g')+r'}{\ohm}'
(ROOT/'datos/parametros.tex').write_text('% Generado por recalcular.py\n'+'\n'.join(
    '\\newcommand{\\'+k+'}{'+v+'}' for k,v in macros_tex.items())+'\n',encoding='utf-8')
(ROOT/'datos/ajustes.json').write_text(json.dumps(FITS,indent=2,ensure_ascii=False),encoding='utf-8')

def exportar_datos_residuos():
    """Entrega las parejas x,y del ajuste, antes del redondeo de presentación."""
    destination=ROOT/'datos/residuos'
    destination.mkdir(exist_ok=True)
    for circuit in ['RC','RL']:
        coordinate='1/f_exp (s)' if circuit=='RC' else 'f_exp (Hz)'
        symbol='C' if circuit=='RC' else 'L'
        lines=[f'# {circuit}: x = {coordinate}; y = X_{symbol} (ohm)',
               '# Dos columnas: x y. Cinco puntos, sin redondear las tablas del informe.']
        for row in RESULTS[circuit]:
            x=1/row['f_Hz'] if circuit=='RC' else row['f_Hz']
            lines.append(format(x,'.17g')+' '+format(row['X'],'.17g'))
        (destination/f'residuos_{circuit}.txt').write_text('\n'.join(lines)+'\n',encoding='utf-8')

    document=r'''\documentclass[10pt,letterpaper]{article}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage[provide=*,spanish]{babel}
\usepackage{lmodern,geometry,amsmath,booktabs,siunitx}
\geometry{margin=2cm}
\setlength{\parindent}{0pt}
\setlength{\parskip}{0.35em}
\sisetup{output-decimal-marker={,},group-digits=false}
\begin{document}
\begin{center}
{\Large\bfseries Datos para el análisis de residuos: RC y RL}
\end{center}
Se conservan los cinco puntos de cada circuito empleados en el informe.
Los archivos de texto contienen dos columnas, $x$ e $y$, y mantienen los valores
previos al redondeo de las tablas del informe. La frecuencia experimental se obtiene
del período medido de la señal del generador.

\section*{Circuito RC}
Para el ajuste lineal se emplean $x=1/f_{\rm exp}$, en segundos,
e $y=X_C$, en ohmios. Archivo: \texttt{residuos\_RC.txt}.
\begin{center}
\small
\begin{tabular}{cccc}
\toprule
$f_{\rm nom}$ (Hz) & $f_{\rm exp}$ (Hz) & $x=1/f_{\rm exp}$ (s) & $y=X_C$ ($\Omega$) \\
\midrule
'''
    for row in RESULTS['RC']:
        values=[row['f_nom'],row['f_Hz'],1/row['f_Hz'],row['X']]
        document+=' & '.join(texnum(value,12) for value in values)+r' \\'+'\n'
    document+=r'''\bottomrule
\end{tabular}
\end{center}

\section*{Circuito RL}
Se emplean $x=f_{\rm exp}$, en hertz, e $y=X_L$, en ohmios.
Archivo: \texttt{residuos\_RL.txt}.
\begin{center}
\small
\begin{tabular}{ccc}
\toprule
$f_{\rm nom}$ (Hz) & $x=f_{\rm exp}$ (Hz) & $y=X_L$ ($\Omega$) \\
\midrule
'''
    for row in RESULTS['RL']:
        values=[row['f_nom'],row['f_Hz'],row['X']]
        document+=' & '.join(texnum(value,12) for value in values)+r' \\'+'\n'
    document+=r'''\bottomrule
\end{tabular}
\end{center}

\section*{Uso con tu cuaderno}
Coloca los dos archivos \texttt{.txt} en la carpeta desde la que ejecutas el cuaderno.
Después de la celda que importa NumPy y define \texttt{graficar\_ajuste}, puedes usar:
\begin{verbatim}
for archivo, titulo in [
    ("residuos_RC.txt", "RC: XC frente a 1/f"),
    ("residuos_RL.txt", "RL: XL frente a f"),
]:
    x, y = np.genfromtxt(archivo).T
    a, b = np.polyfit(x, y, 1)
    residuos = graficar_ajuste(
        x, y, lambda z: a*z + b, titulo
    )
\end{verbatim}
El ajuste $y=ax+b$ utiliza los cinco puntos con intercepto libre.
Los residuos son $r_i=y_i-(ax_i+b)$ y se expresan en ohmios en ambos circuitos.
\end{document}
'''
    (destination/'Datos_Analisis_Residuos_Lab4.tex').write_text(document,encoding='utf-8')

exportar_datos_residuos()
assert len(RAW)==30 and all(len(RESULTS[c])==5 for c in RESULTS)
assert abs(RESULTS['RL'][0]['VG']-1.25)<1e-12 and abs(RESULTS['RL'][0]['VR']-1)<1e-12
assert all(q['a']>0 for q in FITS.values())
print(json.dumps(FITS,indent=2,ensure_ascii=False))
