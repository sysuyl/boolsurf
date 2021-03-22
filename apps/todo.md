# Theory
1. Self intersections: verificare che basta controllare le parita’ del label.
2. Cicli: etichettare celle nel ciclo col poligono entrante, rimuovere gli archi dei poligoni coinvolti nel ciclo, continuare con visita normales sui DAGs rimanenti.

# Implementation
1. Visualizzare booleane.
2. Tirare su tests.
3. Rappresentare il risultato implicitamente come sequenze di puniti.
4. Integrare bezier.

# Performance
1. Farsi a mano triangle-split con singolo segmento (senza chiamare CDT).
2. Provare varie hashmaps.
3. Provare `border_tags` rappresentato in maniera sparse invece che densa, ossia con `hash_map<int, vec3i>` invece che `vector<vec3i>`
4. Parallelismo?    