# Autore: Rodolfo Tolloi

# L'obbiettivo di questo programma è quello di colorare
# la mappa del Canada con 4 colori, assegnando colori diversi
# a regioni adiacenti. Per fare ciò creerò una QUBO 
# composta da due vincoli: il primo è quello di scegliere un
# solo colore per regione, mentre il secondo è di avere
# regioni confinanti con colori diversi. 
# Verrano utilizzate diverse librerie:
# -da dwawe.system verrano importati il sampler e il
# composite per l'embedding
# -dwave.inspector
# -pyqubo, per creare la formulazione matematica del problema
# -pprint, per visualizzare meglio i dati
# -networkx
# -matplotlib

# Importo il compiste per l'embedding e creo il sampler
from dwave.system import EmbeddingComposite, DWaveSampler
sampler = EmbeddingComposite(DWaveSampler(solver={'qpu': True}))

# Analizzatore del d-wave
import dwave.inspector

# Importo da pyqubo le funzioni per creare i problemi
from pyqubo import Array, Binary
from pyqubo import  UserDefinedExpress

# Importo pprint e fisso delle impostazioni adatte
import pprint
pp = pprint.PrettyPrinter(indent=4)

# Importo networkx per creare i grafi
import networkx as nx

# Importo matplotlib per graficare i grafi
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt

# Uso la funzione array di pyqubo per creare dei vettori
# di 4 variabili binarie, le quali rappresentano i colori, per ogni provincia
bc = Array.create('bc', shape=4, vartype='BINARY')   # British Columbia
ab = Array.create('ab', shape=4, vartype='BINARY')   # Alberta
sk = Array.create('sk', shape=4, vartype='BINARY')   # Saskatchewan
mb = Array.create('mb', shape=4, vartype='BINARY')   # Manitoba
on = Array.create('on', shape=4, vartype='BINARY')   # Ontario
qc = Array.create('qc', shape=4, vartype='BINARY')   # Quebec
nl = Array.create('nl', shape=4, vartype='BINARY')   # Newfoundland and Labrador
nb = Array.create('nb', shape=4, vartype='BINARY')   # New Brunswick
pe = Array.create('pe', shape=4, vartype='BINARY')   # Prince Edward Island
ns = Array.create('ns', shape=4, vartype='BINARY')   # Nova Scotia
yt = Array.create('yt', shape=4, vartype='BINARY')   # Yukon
nt = Array.create('nt', shape=4, vartype='BINARY')   # Northwest Territories
nu = Array.create('nu', shape=4, vartype='BINARY')   # Nunavut  

# Creo un vettore con all'interno tutte le provincie 
provinces = [bc, ab, sk, mb, on, qc, nl, nb, pe, ns, yt, nt, nu]


# Creo un vettore con le coppie di provincie vicine
neighbours = [(bc, ab),
              (bc, nt),
              (bc, yt),
              (ab, sk),
              (ab, nt),
              (sk, mb),
              (sk, nt),
              (mb, on),
              (mb, nu),
              (on, qc),
              (qc, nb),
              (qc, nl),
              (nb, ns),
              (yt, nt),
              (nt, nu)]

# Creo una funzione che, dato un vettore di variabili, crea
# l'equazione che rappresenta il vincolo di scelta di un solo colore
class one_color(UserDefinedExpress):
 def __init__(self, a):
    express = (a[0]+a[1]+a[2]+a[3]-1)**2
    super().__init__(express)

# Creo la prima parte dell'Hamiltoniana unendo i vincoli di tutte 
# le 13 provincie
H1=one_color(provinces[len(provinces)-1])
for i in range(len(provinces)-1):
    H1=H1+one_color(provinces[i])

# Creo la funzione che, data una coppia di vettori, crea la funzione
# che penalizza il fatto che questi siano dello stesso colore
class colori_diversi(UserDefinedExpress):
 def __init__(self, a, b):
    express = a[0]*b[0]+a[1]*b[1]+a[2]*b[2]+a[3]*b[3]
    super().__init__(express)

# Creo la seconda parte dell'Hamiltoniana creando e unendo 
# tutti i vincoli per le coppie di provincie vicine
H2=colori_diversi(neighbours[len(neighbours)-1][0],neighbours[len(neighbours)-1][1])
for i in range(len(neighbours)-1):
    H2=H2+colori_diversi(neighbours[i][0],neighbours[i][1])

# Creo l'Hamiltoniana e la faccio diventare una BQM
model = H1+H2
bqm = model.compile().to_bqm()


# Creo il set dei risultati utilizzando il sampler definito prima
chainstrenght= None
sampleset = sampler.sample(bqm, num_reads = 10,
     chain_strength=chainstrenght, label="Problema della mappa")

#risultati
#print(sampleset)

# Uso l'inspector
dwave.inspector.show(sampleset)

# Seleziono il risultato migliore, ossia quello a energia minore
risultato = sampleset.first.sample

#print(risultato)
#print(risultato[0])

# Creo un grafo vuoto
G = nx.Graph()

# Creo un vettore con i nodi e uno con gli archi
nodes=[ 'bc',
        'ab',
        'sk',
        'mb',
        'on',
        'qc',
        'nl',
        'nb',
        'pe',
        'ns',
        'yt',
        'nt',
        'nu']

edges =     [('bc', 'ab'),
              ('bc', 'nt'),
              ('bc', 'yt'),
              ('ab', 'sk'),
              ('ab', 'nt'),
              ('sk', 'mb'),
              ('sk', 'nt'),
              ('mb', 'on'),
              ('mb', 'nu'),
              ('on', 'qc'),
              ('qc', 'nb'),
              ('qc', 'nl'),
              ('nb','ns'),
              ('yt', 'nt'),
              ('nt', 'nu')]

# Creo il grafo che rappresenta la mappa del Canada              
G.add_nodes_from(nodes)
G.add_edges_from(edges)

# Creo un vettore i cui elementi sono il nome della provincia
# e il colore assegnato, sotto forma di numero
color_labels = [k for k, v in risultato.items() if v == 1]

# Divido l'elemento del vettore in nome e colore
for i in range(len(color_labels)):
        name = color_labels[i][:2]
        color = color_labels[i][3]
        G.nodes[name]["color"] = color

# Creo un vettore di dimensioni 13 con ogni colore associato alla provincia,
# sempre in forma di numero
color_map = [color for name, color in G.nodes(data="color")]

# Sostituisco il numero che indica il colore con la stringa
# intepretabile da networkx 
for i in range(13):
    color_map[i]=color_map[i].replace('0','red')
    color_map[i]=color_map[i].replace('1','blue')
    color_map[i]=color_map[i].replace('2','green')
    color_map[i]=color_map[i].replace('3','violet')
    
# Creo un vettore con le posizioni dei nodi nel grafico
node_positions = {"bc": (0, 1),
                  "ab": (2, 1),
                  "sk": (4, 1),
                  "mb": (6, 1),
                  "on": (8, 1),
                  "qc": (10, 1),
                  "nb": (10, 0),
                  "ns": (12, 0),
                  "pe": (12, 1),
                  "nl": (12, 2),
                  "yt": (0, 3),
                  "nt": (2, 3),
                  "nu": (6, 3)}

# Creo il grafico del grafo
nx.draw_networkx(G, pos=node_positions, with_labels=True,
                   node_color=color_map, font_color="w", node_size=400)

# Salvo il grafo in un file
filename = "mappa-colorata.png"
plt.savefig(filename)
print("The graph is saved in '{}'.".format(filename))

