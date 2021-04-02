# Algorithm
- [x] Self intersections: verificare che basta controllare le parita’ dei label.
- [x] Cicli: etichettare celle nel ciclo, rimuovere gli archi dei poligoni coinvolti nel ciclo, continuare con visita standard sui DAGs rimanenti.
- [x] Trovare quale fra le celle ambiente candidate e' quella giusta.

# Implementation
- [x] Visualizzare booleane.
- [x] Tirare su tests.
- [x] Rappresentare il risultato implicitamente come sequenze di punti: inserire intersezioni nei poligoni e taggare i lati per visualizzare bene le booleane.
- [ ] Integrare bezier.

# Performance
- [x] Farsi a mano triangle-split con singolo segmento (senza chiamare CDT).
- [ ] Provare varie hashmaps.
- [x] Riordinare triangoli per flood fill.
- [ ] Provare border_tags rappresentato in maniera sparsa invece che densa, ossia con `hash_map<int, vec3i>` invece che `vector<vec3i>`

# Misc
- [x] Script svg parser con suddivisioni da cli
- [x] Chiamare svg_parser.py da c++ per convertire .svg come cli 
- [ ] Triangolatore di Cino 
- [ ] Fare immagini belle
- [ ] Runnare test 
  
# Giacomo
- [x] Fix windows build (freetype)
- [ ] Fix windows rendering

# Marzia
- [x] Rimuovere segments da polygon (lasciare solo edges ovunque)
- [x] Fix compilazione in release windows (?)
