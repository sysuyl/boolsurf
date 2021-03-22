# Grafo delle celle
Un nodo per ogni cella. Ogni cella e' connessa a tutte le celle a cui e' adiacente sulla superficie.
L'arco e' etichettato con le seguenti informazioni:
- Cella a cui arriva l'arco
- Poligono che viene attversato andando dalla cella a quella adiacente
- Orientamento del poligono...

# Cella ambiente
Sulle superfici non e' detto che ci sia una cella ambiente, ossia che e' interna a nessun poligono.
Quando siamo sicuri che ne esiste una (vedi "Cicli"), definiamo la cella ambiente come la cella 
senza archi entranti. Possono esserci anche piu' di una cella ambiente, ma per funzionare l'algoritmo
basta che ne troviamo una. Il labeling della cella ambiente viene settato a "0" e inizia la visita
da li'. In realta' neanche ci serve farlo, grazie al fix di [Label negativi] possiamo partire da qualsiasi
cella e aggiustare i labeling in un secondo pass. Infatti partire dalla cella sbagliata, fa sbagliare
i labaling a meno di una costante. 

# Percorsi multipli
Nello spazio euclideo i grafi sono "conservativi": date due celle ogni percorso da una all'altra da' lo stesso risutltato. In altre topologie, dipende dal percorso. Vedi esempio percorsi.json

In caso di percorsi multipli, prendiamo il massimo dei labeling che "arrivano"

Una conseguenza dei percorsi multipli e' che possono esserci anche cicli.

# Cicli
Troviamo i cicli. Troviamo i poligoni che fanno parte dei cicli. Settiamo il labeling iniziale
delle celle che fanno parte dei cicli. Facciamo partire la visita, partendo dalle celle che fanno
parte dei cicli, ignorando tutti gli archi nell'adiacneza relativi ai poligoni che fanno parte dei
cicli. In questo caso non cerchiamo neanche le celle ambiente.


# Self-interesections
Nel caso self-intersections si possono generare labeling > 1. In quei casi facciamo % 2.
Praticmente facciamo even-odd rule.

# Label negativi
Nel caso in cui la visita produca label negativi, troviamo il minimo offset per label che faccia in 
modo che, se aggiunto a tutti, non rimanga nessun label negativo.


# Pesudocodice
cicli := calcola_cicli()