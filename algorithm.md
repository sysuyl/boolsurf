# Immagini
1. Modello con poligoni sopra e intersezioni visualizzate.
2. Modello triangolato lungo i bordi dei poligoni (si vedono i mesh-edges)
3. Celle colorate di colori diversi (flood-filll).
4. Aggiungiamo labeling in sovraimpressione alle celle come testo. Celle con lo stesso labeling diventano dello stesso colore.
5. Mostriamo risultato dell’operazione booleana, disegnando i control points dei nuovi poligoni risultanti.

# 1 [slice_mesh]
- Add polygon points as mesh vertices (snapping to original vertices if too close).
- Creiamo hashgrid: ad ogni faccia della mesh associamo una lista di polilinee (?).
  - (Immagine piccolina per spiegare cosa sono le polilinee)
- Aggiungiamo intersezioni all'hashgrid
  
# 2 [slice_mesh]
- Per ogni faccia dell'hashgrid
  - Calcoliamo i nodi di input per CDT
  - Calcoliamo i vincoli di edge per CDT
  - Triangoliamo con CDT che restituisce triangoli e le loro adianceze (locali)
  - Usiamo i triangoli locali per aggiornare i triangoli globali 
  - Usiamo le adiacenze locali per aggiornare le adiacenze globali
  - Aggiustiamo le adiancenze dei triangoli adiancenti a quelli triangolati [update_face_adjacencies]
  - Calcoliamo border-tags
    - Creiamo hashmap che mappa gli edge che embeddano i poligoni al poligono di appartenenza 
    - Per ogni triangolo triangolato e per ogni edge di quel triangolo
      - Se l'edge sta nell'hashmap, tagghiamo quel lato del triangolo come interno al poligono (che è il valore associato all'edge nell'hashmap)
      - Se l'edge al contratio sta nell'hashmap, tagghiamo quel lato del triangolo come esterno al poligono (che è il valore associato all'edge nell'hashmap)      -  
  - 

# 3 [flood-fill]
- Il flood-fill trova le celle e genera il grafo delle adiacenze tra esse. In questo grafo c'è un nodo per ogni cella. Ogni cella è connessa a tutte le celle a cui è adiacente sulla superficie.
L’arco è etichettato con le seguenti informazioni:
- Poligono che viene attraversato andando dalla cella a quella adiacente
- Un bit che rappresenta l'orientamento del poligono attraversato (ovvero se sto entrando o uscendo da esso)
    
- Inizializziamo lo stack con l'ultima faccia della mesh
- Finchè lo stack è pieno:
  - Estaiamo l'ultima faccia dallo stack
  - Se la faccia è già stata taggata come appartenente ad una cella, allora skippiamo
  - Aggiungiamo una nuova cella al risultato
  - Ci espandiamo tramite il flood-fill fino a che non troviamo un bordo taggato (che rappresenta una possibile nuova cella)

    neighbor_cell = tag della faccia adiacente all'edge corrente
    if (neighbor_cell >= 0) //Se la cella è stata già visitata
        if (neighbor_cell == cell_id)
            if (crossing border)
                connetti cella a se stessa tramite l'arco di valore +p
            else
                continue //ho già visitato durante il flood-fill di questa cella
        else
            colleghiamo cella a neighbor
    else //se non ho mai visitato la faccia adiacente
        if (crossing border)
            aggiungo neighbor al cell_stack //mi espanderò su questa faccia dopo aver completato il flood-fill di questa cella
        else
            aggiungo neighbor al face stack //mi espanderò su questa faccia

# 4 [compute_cell_labels] and [propagate_cell_labels]
- Calcoliamo i cicli presenti nel grafo delle adiacenze tra celle.
- Per ogni arco coinvolto in un ciclo tagghiamo la cella destinazione dell'arco come interna al poligono identificato dall'arco. 
- Eliminiamo tutti gli archi taggati con i poligoni presenti nei cicli. 
- Propaghiamo i label nel grafo di adiacenza delle celle. Se nel grafo non sono presenti cicli, partiamo da una qualsiasi delle celle. Se ci sono cicli, iniziamo la visita a partire dai noi coinvolti nei cicli. 
  - La propagazione è analoga a mesh arrangements e segue solo gli archi positivi
  - Se un nodo viene raggiunto durante la visita con labelling diversi allora prendiamo il labelling massimo per ogni componente
- Se sono rimasti nodi non visitati facciamo ripartire la visita dai nodi che hanno almeno un vicino non visitato e seguendo solo gli archi negativi
- Al termine di queste operazioni alcuni labelling potrebbero essere negativi. Calcoliamo un labelling offset tale che se sommato a tutti i labelling nessuno risulta più negativo
- Nel caso in cui la label sia > 1 applichiamo la even-odd rule
  - Nelle self-intersections posso entrare in un poligono più volte senza esserne prima uscito)

# 5 [compute_bool_operation]
- Le operazioni booleane tra due poligoni vengono effettuate con una selezione delle celle che li compongono in base ai labelling
- 
        



