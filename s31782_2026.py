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
