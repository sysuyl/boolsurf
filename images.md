# Paper Images
- Teaser formato da 3 immagini di cui
  - una con i poligoni tracciati, evidenziando i punti di controllo e le intersezioni
  - una visualizzazione delle celle calcolate e di un'operazione booleana
  - una con i poligoni risultanti dalle booleane 

- Sequenza di immagini che spiega l'algoritmo (vedi anche algorithm.md)
  - Modello con i poligoni iniziali e con intersezioni disegnate
  - Triangolazione (visualizzare bene gli edge embeddati)
  - Identificazione delle celle con colori differenti
  - Celle con la stessa label sono colorate allo stesso modo e con i label in sovraimpressione
  - Risultato finale della booleana con annessi control-point e brodi 
  
- Spiegazione della triangolazione e della creazione dei vincoli
  - Immagine in cui facciamo vedere cosa accade all'intendo di un triangolo (+ i suoi vicini?), rappresentando quali sono gli edgeaggiunti alla triangolazione con CDT

- Spiegazione del flood-fill
  - Sequenza di immagini che fanno vedere come si espande il flood-fill (e in contemporanea come si crea il grafo delle adiacenze)

- Labelling (visita del frago delle adiancenze)
  - Immagine per far capire la gestione dei cicli
  - Immagini per far capire come avviene la visita del grafo e per spiegare i casi particolari (funzionamento normale, self-loops, cicli e altri casi che mo non mi ricordo) 

- Risultati, prendere modelli belli da thingi10k o da altri paper (Crane!) e disegnare a mano diversi esempi tipo
  - Tante intersezioni tra poligoni
  - Numero alto di suddivisioni delle curve
  - Poligoni che passano su tanti triangoli
  - Poligoni che passano su un triangolo solo

- Test su Thingi10k (da mostrare a parte nei supplemental - 4321 mesh totali)
  - Qui va capito se mostrare il risultato di una booleana come area colorata o solamente tramite il bordo 


