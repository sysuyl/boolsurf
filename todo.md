# Theory
- [x] Self intersections: verificare che basta controllare le parita’ dei label.
- [ ] Cicli: etichettare celle nel ciclo, rimuovere gli archi dei poligoni coinvolti nel ciclo, continuare con visita standard sui DAGs rimanenti.
- [ ] Trovare quale fra le celle ambiente candidate e' quella giusta.

# Implementation
- [ ] Visualizzare booleane.
- [ ] Tirare su tests.
- [ ] Rappresentare il risultato implicitamente come sequenze di puniti.
- [ ] Integrare bezier.

# Performance
- [ ] Farsi a mano triangle-split con singolo segmento (senza chiamare CDT).
- [ ] Provare varie hashmaps.
- [ ] Provare border_tags rappresentato in maniera sparsa invece che densa, ossia con `hash_map<int, vec3i>` invece che `vector<vec3i>`
- [ ] Parallelismo?

# Giacomo
- Fix windows build
- Fix windows rendering
