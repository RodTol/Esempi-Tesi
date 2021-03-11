# Autore: Rodolfo Tolloi

# Con questo programma voglio creare un codice che, dato un grafo
# qualsiasi, trovi il suo maximum indipendet set, ossia il numero
# massimo di nodi non connessi all'interno del grafo.
# Per fare ciò, andrà per prima cosa creata una QUBO che abbia
# come obbiettivo quello di individuare con la minima energia il set
# indipendete più grande.
# Succesivamente cercherè i risultati usando uno dei sample del D-wawe,
# ripetendo l'annealing più volte e scegliendo come risultato migliore
# quello con l'energia più bassa
# Infine riconverirò il risultato migliore da una list a un
# grafo che proietterò sopra a quello orginale. 
# Per creare questo procedimento ho modificato la funzione
# dnx.maximum_independent_set dal pacchetto dwawe_networkx.
#
# Librerie:
# - networkx per creare visualizzare e modificare i grafi
# - il pacchetto matplotlib.pyplot per graficare
# - dwawe.inspector 
# - dwawe.system.sampler
# - dwave.system.composites

# Importo il pacchetto per modificare i grafi
import networkx as nx

# Importo la libreria per graficare i grafi
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt

# Importo l'analizzatore
import dwave.inspector

# Selezioni il sampler che useremo. Voglio un metodo quantistico,
# quindi uso il DwaweSampler
from dwave.system.samplers import DWaveSampler

# Utilizzo un composite per fare l'embedding
from dwave.system.composites import EmbeddingComposite

# Preparo il sampler usando l'embedding scelto
sampler = EmbeddingComposite(DWaveSampler(solver='Advantage_system1.1'))

# Creo un grafo vuoto
G = nx.Graph()

# Prendo un grafo da un file
grafo = open("JOHNSON8-2-4.txt", "rb")
#grafo = open("/workspace/Tesi/Esempi-commentati/Pipelines-Antennas/chesapeake.txt", "rb")
G = nx.read_edgelist(grafo)

#Creo la funzione per formulare il QUBO 
def massimo_set_indipendente_qubo(G, weight=None, lagrange=2.0):
# un QUBO vuoto per un grafo vuoto
    if not G:
        return {}

    # Definiamo il set indipendente più largo come S.
    # Per ogni nodo n nel grafo, assegniamo una variabile binaria, v_n, che
    # sarà v_n=1 se n è in S, altrimenti v_n=0.
    # Chiamiamo la matrice del problema QUBO, Q.
    # Sulla diagonale di questa poniamo un bias linear, uguale al negativo del peso
    # di quel nodo. In questo modo, ogni nodo sarà indirizzato a essere in S.
    # I pesi sono normalizzati per essere al massimo uguali a 1. Pesi negativi
    # sono considerati nulli.
    # I termini non diagonali sono invece tutti uguali a 2, in maniera tale che
    # se entrambi i nodi appartengono a S, l'energia totale sarà aumentata di 2.
    cost = dict(G.nodes(data=weight, default=1))
    scale = max(cost.values())
    Q = {(node, node): min(-cost[node] / scale, 0.0) for node in G}
    Q.update({edge: lagrange for edge in G.edges})

    return Q

Q = massimo_set_indipendente_qubo(G, weight = None, lagrange=2.0)


#Faccio l'annealing per 100 volte e creo dunque il mio sampleset
response = sampler.sample_qubo(Q,chain_strength=5,
     num_reads=50, label='Problema delle antenne')

# Seleziono il risultato con l'energi più bassa e lo visualizzo
sample = next(iter(response))
# print(sample)

# I nodi che hanno spin up sono quelli presenti nel set minimo
# mentre quelli con spin nullo no. Creo quindi la mia lista di nodi
S = [node for node in sample if sample[node] > 0]

# Uso l'inspector sul mio sampleset
dwave.inspector.show(response)

# Visualizzo il risultato ottenuto
print('Maximum independent set size found is', len(S))
print(S)

# Vado a creare il grafo con i nodi trovati dalla soluzione
k = G.subgraph(S)
notS = list(set(G.nodes()) - set(S))
othersubgraph = G.subgraph(notS)
pos = nx.spring_layout(G)
plt.figure()

# Salvo e stampo una immagine con il grafo originale
original_name = "antenna_plot_original.png"
nx.draw_networkx(G, pos=pos, with_labels=True)
plt.savefig(original_name, bbox_inches='tight')

# Salvo la soluzione e la stampo
# NOTA: i nodi in rosso sono quelli che fanno parte del mio maximum set,
# quelli in blu no
solution_name = "antenna_plot_solution.png"
nx.draw_networkx(k, pos=pos, with_labels=True, node_color='r', font_color='k')
nx.draw_networkx(othersubgraph, pos=pos, with_labels=True, node_color='b', font_color='w')
plt.savefig(solution_name, bbox_inches='tight')

print("I plot sono salvati in {} e {}".format(original_name, solution_name))
