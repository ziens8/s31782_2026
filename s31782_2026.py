import random
import csv


# ---------- Funkcje wymagane przez specyfikację ----------

def generate_sequence(length: int) -> str:
    """Zwraca losową sekwencję DNA o zadanej długości.

    Korzysta z random.choices, które zwraca listę k losowo wybranych
    elementów z populacji nukleotydów {A, C, G, T} z rozkładem jednostajnym.
    """
    nucleotides = ['A', 'C', 'G', 'T']
    # ''.join scala listę pojedynczych znaków w jeden string
    return ''.join(random.choices(nucleotides, k=length))


def calculate_stats(sequence: str) -> dict:
    """Zwraca słownik ze statystykami sekwencji.

    Klucze:
        "A", "C", "G", "T" - procentowy udział danego nukleotydu (float, %),
        "GC"               - zawartość GC w procentach (float, %).
    """
    n = len(sequence)
    if n == 0:
        # zabezpieczenie przed dzieleniem przez zero (teoretycznie nieosiągalne
        # bo walidator nie pozwoli na długość 0, ale dobre obyczaje)
        return {"A": 0.0, "C": 0.0, "G": 0.0, "T": 0.0, "GC": 0.0}

    stats = {}
    for nt in ("A", "C", "G", "T"):
        # str.count(x) zwraca liczbę wystąpień podciągu x w stringu
        stats[nt] = sequence.count(nt) / n * 100

    # GC-content to po prostu suma udziałów G i C
    stats["GC"] = stats["G"] + stats["C"]
    return stats


def insert_name(sequence: str, name: str) -> str:
    """Wstawia imię w losową pozycję sekwencji. Imię zapisane małymi literami."""


    if not name:
        return sequence  # brak imienia -> nic nie wstawiamy

    # randint(a, b) zwraca liczbę z przedziału [a, b] - dzięki len(sequence)
    # jako górnej granicy dopuszczamy wstawienie na samym początku lub końcu
    pos = random.randint(0, len(sequence))
    return sequence[:pos] + name.lower() + sequence[pos:]


def format_fasta(seq_id: str, description: str,
                 sequence: str, line_width: int = 80) -> str:
    """Zwraca sformatowany rekord FASTA jako pojedynczy string.

    Pierwsza linia: 'ID opis' (jeśli opis pusty, spacja jest usuwana przez rstrip).
    Kolejne linie: sekwencja łamana co line_width znaków.
    """
    header = f">{seq_id} {description}".rstrip()

    # list comprehension dzieli sekwencję na kawałki długości line_width:
    # dla 'ABCDEFGHIJ' i line_width=4 dostaniemy ['ABCD', 'EFGH', 'IJ']
    lines = [sequence[i:i + line_width]
             for i in range(0, len(sequence), line_width)]

    return header + "\n" + "\n".join(lines)


def validate_positive_int(prompt: str,
                          min_val: int = 1,
                          max_val: int = 100_000) -> int:
    """Pobiera od użytkownika liczbę całkowitą z zakresu [min_val, max_val].

    Przy błędnym wejściu (nie-liczba, poza zakresem) powtarza pytanie zamiast
    rzucać wyjątkiem. Pętla działa, aż otrzyma poprawną wartość.
    """
    while True:
        raw = input(prompt)
        try:
            value = int(raw)
        except ValueError:
            # int() rzuca ValueError gdy wejście nie da się skonwertować
            print(f"Błąd: wartość musi być liczbą całkowitą "
                  f"z zakresu [{min_val}, {max_val}].")
            continue

        if value < min_val or value > max_val:
            print(f"Błąd: wartość musi być liczbą całkowitą "
                  f"z zakresu [{min_val}, {max_val}].")
            continue

        return value

# ---------- Funkcja pomocnicza ----------

def validate_id(prompt: str) -> str:
    """Pobiera ID sekwencji - nie może być puste ani zawierać białych znaków."""
    while True:
        seq_id = input(prompt).strip()
        if not seq_id:
            print("Błąd: ID nie może być puste.")
            continue
        # split() bez argumentów dzieli po dowolnych białych znakach;
        # jeśli wynik ma więcej niż 1 element, znaczy że w ID była spacja/tab
        if len(seq_id.split()) > 1:
            print("Błąd: ID nie może zawierać białych znaków.")
            continue
        return seq_id
# ---------- Dodatki (4 sztuki, wymagane min. 4) ----------

def reverse_complement(sequence: str) -> tuple[str, str]:
    """Zwraca krotkę (komplementarna, odwrotnie_komplementarna).

    Parowanie zasad: A<->T, C<->G.
    Nić komplementarna - czytana w tym samym kierunku co oryginał.
    Nić odwrotnie komplementarna - komplementarna odczytana 5'->3', czyli
    odwrócona (w biologii to faktyczna druga nić DNA).
    """
    pairing = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    complement = ''.join(pairing[nt] for nt in sequence)
    reverse_comp = complement[::-1]  # [::-1] to idiom Pythonowy na odwrócenie
    return complement, reverse_comp


def transcribe(sequence: str) -> str:
    """Transkrypcja in silico: DNA -> mRNA przez zamianę T na U."""
    return sequence.replace('T', 'U')


def find_motif(sequence: str, motif: str) -> list[int]:
    """Zwraca listę pozycji (indeksowanie od 1) wystąpień motywu.

    Wyszukuje też nakładające się wystąpienia (w przeciwieństwie do str.find
    w pętli, która by je pomijała).
    """
    positions = []
    if not motif:
        return positions

    motif = motif.upper()
    m = len(motif)
    # iterujemy po wszystkich możliwych pozycjach startowych okna
    for i in range(len(sequence) - m + 1):
        if sequence[i:i + m] == motif:
            positions.append(i + 1)  # +1 bo biolodzy liczą od 1, nie od 0
    return positions


def sliding_window_gc(sequence: str, window: int) -> list[tuple[int, float]]:
    """GC-content w oknie przesuwnym o szerokości `window`, krok = 1.

    Zwraca listę krotek (pozycja_startu_od_1, gc_content_w_procentach).
    """
    results = []
    if window <= 0 or window > len(sequence):
        return results

    for i in range(len(sequence) - window + 1):
        win = sequence[i:i + window]
        gc = (win.count('G') + win.count('C')) / window * 100
        results.append((i + 1, round(gc, 2)))
    return results


def save_sliding_window_csv(data: list[tuple[int, float]], filename: str) -> None:
    """Zapisuje dane sliding window do pliku CSV (kolumny: pozycja, gc_content)."""
    with open(filename, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["pozycja_startu", "gc_content"])
        writer.writerows(data)
