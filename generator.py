import yaml
import sympy as sp
from sympy import simplify
from collections import defaultdict
from sympy import sympify


def load_yaml(file_path):
    """
    Lädt eine YAML-Datei und gibt die enthaltenen Daten als Dictionary zurück.
    :param file_path: Pfad zur YAML-Datei.
    :return: Dictionary mit den geladenen Daten.
    """
    with open(file_path, "r") as file:
        data = yaml.safe_load(file)
    return data

def get_slow_educts(reaction, slow_species):
    """
    Extrahiert die Koeffizienten der Edukte, die zu den langsamen Spezies gehören.
    :param reaction: Ein Dictionary, das eine Reaktion beschreibt.
    :param slow_species: Eine Liste der langsamen Spezies.
    :return: Dictionary mit den Koeffizienten der langsamen Edukte.
    """
    educts = reaction.get("educts", {})
    slow_educts = {species: coeff for species, coeff in educts.items() if species in slow_species}
    return slow_educts

def init_slow_symbolic_variables(slow_species):
    """
    Erstellt symbolische Variablen für jede langsame Spezies.
    :param slow_species: Liste der langsamen Spezies.
    :return: Dictionary mit den symbolischen Variablen.
    """
    return {species: sp.Symbol(f'v_{{{species}}}') for species in slow_species}

def init_slow_symbolic_derivatives(slow_species):
    """
    Erstellt symbolische Variablen für die Ableitungen nach den langsamen Spezies.
    :param slow_species: Liste der langsamen Spezies.
    :return: Dictionary mit den symbolischen Ableitungsvariablen.
    """
    return {species: sp.Symbol(f'df_{{{species}}}') for species in slow_species}

def multiply_slow_species(reaction, slow_symbolic_variables):
    """
    Multipliziert die symbolischen langsamen Speziesvariablen mit den Exponenten aus den Edukten der Reaktion.
    :param reaction: Ein Dictionary, das eine Reaktion beschreibt.
    :param slow_symbolic_variables: Dictionary der symbolischen Variablen für langsame Spezies.
    :return: Symbolischer Ausdruck der multiplizierten Variablen.
    """
    slow_educts = get_slow_educts(reaction, slow_symbolic_variables.keys())
    return sp.Mul(*(slow_symbolic_variables[species]**exp for species, exp in slow_educts.items())) if slow_educts else 1

def sum_slow_educt_coefficients(reaction, slow_species):
    """
    Summiert die Koeffizienten der langsamen Edukte einer Reaktion.
    :param reaction: Ein Dictionary, das eine Reaktion beschreibt.
    :param slow_species: Liste der langsamen Spezies.
    :return: Summe der Koeffizienten der langsamen Edukte.
    """
    slow_educts = get_slow_educts(reaction, slow_species)
    return sum(slow_educts.values())

def compute_scaled_rate(reaction, natnum, slow_species):
    """
    Berechnet die skalierte Rate der Reaktion: rate * (N ** (scale - 1 + sum_slow_educt_coefficients(reaction))).
    :param reaction: Ein Dictionary, das eine Reaktion beschreibt.
    :param natnum: Symbolische Variable für N.
    :param slow_species: Liste der langsamen Spezies.
    :return: Symbolischer Ausdruck der skalierten Rate.
    """
    rate = sp.Symbol(reaction["rate"])
    scale = sp.Symbol(reaction["scale"])
    f_reaction = sum_slow_educt_coefficients(reaction, slow_species)
    g_reaction = scale - 1 + f_reaction
    return rate * (natnum ** g_reaction)

def compute_species_difference(reaction, species):
    """
    Berechnet die Differenz aus Produkten und Edukten für eine gegebene langsame Spezies.
    :param reaction: Ein Dictionary, das eine Reaktion beschreibt.
    :param species: Die zu betrachtende langsame Spezies.
    :return: Die Differenz aus Produkt- und Edukt-Koeffizienten.
    """
    educts = reaction.get("educts", {}).get(species, 0)
    products = reaction.get("products", {}).get(species, 0)
    return products - educts

def is_slow_only_reaction(reaction, fast_species): 
    """
    Prüft, ob eine Reaktion ausschließlich langsame Spezies konsumiert und produziert.
    :param reaction: Ein Dictionary, das eine Reaktion beschreibt.
    :param fast_species: Liste der schnellen Spezies.
    :return: True, wenn die Reaktion keine schnellen Spezies verbraucht oder produziert, sonst False.
    """
    educts = reaction.get("educts", {}).keys()
    products = reaction.get("products", {}).keys()
    return all(species not in fast_species for species in educts) and all(species not in fast_species for species in products)

def produces_fast_species(reaction, fast_species, all_fast_species):
    """
    Prüft, ob die gegebene schnelle Spezies in den Produkten der Reaktion vorkommt
    und keine schnelle Spezies aus der Liste aller schnellen Spezies in den Edukten enthalten ist.

    :param reaction: Ein Dictionary, das eine Reaktion beschreibt.
    :param fast_species: Die spezifische schnelle Spezies, die geprüft wird.
    :param all_fast_species: Liste aller schnellen Spezies.
    :return: True, wenn die schnelle Spezies produziert, aber keine andere schnelle Spezies verbraucht wird, sonst False.
    """
    educts = reaction.get("educts", {}).keys()
    products = reaction.get("products", {}).keys()

    # Prüft, ob die spezifische fast_species produziert wird UND keine der all_fast_species verbraucht wird
    return fast_species in products and all(species not in all_fast_species for species in educts)

def consumes_fast_species(reaction, fast_species):
    """
    Prüft, ob eine gegebene schnelle Spezies in den Edukten der Reaktion vorkommt.
    :param reaction: Ein Dictionary, das eine Reaktion beschreibt.
    :param fast_species: Liste der schnellen Spezies.
    :return: True, wenn eine schnelle Spezies in den Edukten enthalten ist, sonst False.
    """
    educts = reaction.get("educts", {}).keys()
    return any(species in fast_species for species in educts)

def consumes_fast_species_without_producing(reaction, fast_species):
    """
    Prüft, ob eine schnelle Spezies in den Edukten der Reaktion vorkommt und nicht in den Produkten enthalten ist.
    :param reaction: Ein Dictionary, das eine Reaktion beschreibt.
    :param fast_species: Liste der schnellen Spezies.
    :return: True, wenn eine schnelle Spezies konsumiert aber nicht produziert wird, sonst False.
    """
    educts = reaction.get("educts", {}).keys()
    products = reaction.get("products", {}).keys()
    return any(species in fast_species for species in educts) and all(species not in products for species in fast_species)

def consumes_and_produces_fast_species(reaction, consumed_species, produced_species):
    """
    Prüft, ob eine schnelle Spezies in den Edukten der Reaktion vorkommt und eine andere schnelle Spezies in den Produkten enthalten ist.
    :param reaction: Ein Dictionary, das eine Reaktion beschreibt.
    :param consumed_species: Die schnelle Spezies, die verbraucht werden soll.
    :param produced_species: Die schnelle Spezies, die erzeugt werden soll.
    :return: True, wenn die erste Spezies konsumiert und die zweite Spezies produziert wird, sonst False.
    """
    educts = reaction.get("educts", {}).keys()
    products = reaction.get("products", {}).keys()
    return consumed_species in educts and produced_species in products

def create_fast_species_matrix(reactions, fast_species, slow_symbolic_variables, natnum, slow_species):
    """
    Erstellt eine quadratische Matrix M mit Indizes S, T für schnelle Spezies.

    M[S][S] ist die Summe der Terme für Reaktionen, die S verbrauchen, aber keine andere schnelle Spezies erzeugen.
    M[S][T] (S ≠ T) ist die negative Summe der Terme für Reaktionen, die S verbrauchen und T erzeugen.

    :param reactions: Liste der Reaktions-Dictionaries.
    :param fast_species: Liste der schnellen Spezies.
    :param slow_symbolic_variables: Symbolische Variablen für langsame Spezies.
    :param natnum: Symbolische Variable für N.
    :param slow_species: Liste der langsamen Spezies.
    :return: Dictionary mit symbolischen Matrixeinträgen M[S][T].
    """
    # Initialisiere die Matrix als verschachteltes Dictionary mit 0-Einträgen
    M = {S: {T: sp.S(0) for T in fast_species} for S in fast_species}

    # Berechne die Diagonalelemente M[S][S]
    for S in fast_species:
        for reaction in reactions:
            if consumes_fast_species_without_producing(reaction, S):
                term = (
                    multiply_slow_species(reaction, slow_symbolic_variables) *
                    compute_scaled_rate(reaction, natnum, slow_species)
                )
                M[S][S] += term  # Summiere auf die Diagonale

    # Berechne die Nicht-Diagonal-Elemente M[S][T] (S ≠ T)
    for S in fast_species:
        for T in fast_species:
            if S != T:
                for reaction in reactions:
                    if consumes_and_produces_fast_species(reaction, S, T):
                        term = (
                            multiply_slow_species(reaction, slow_symbolic_variables) *
                            compute_scaled_rate(reaction, natnum, slow_species)
                        )
                        M[S][T] -= term  # Negative Summe für M[S][T]

    return M

def create_modified_fast_species_matrix(reactions, fast_species, slow_symbolic_variables, natnum, slow_species, S, S_prime):
    """
    Erstellt eine modifizierte quadratische Matrix M mit Indizes S', T für schnelle Spezies.

    Die Matrix basiert auf der von create_fast_species_matrix erzeugten Matrix,
    allerdings werden die Einträge M[T][S_prime] wie folgt ersetzt:

        sum_term = sum(
            multiply_slow_species(reaction, slow_symbolic_variables) *
            compute_scaled_rate(reaction, natnum, slow_species) *
            compute_species_difference(reaction, S)
            for reaction in reactions if consumes_fast_species(reaction, T)
        )

    :param reactions: Liste der Reaktions-Dictionaries.
    :param fast_species: Liste der schnellen Spezies.
    :param slow_symbolic_variables: Symbolische Variablen für langsame Spezies.
    :param natnum: Symbolische Variable für N.
    :param slow_species: Liste der langsamen Spezies.
    :param S: Eine einzelne langsame Spezies.
    :param S_prime: Eine einzelne schnelle Spezies.
    :return: Modifizierte symbolische Matrix M[S'][T].
    """
    # Erzeuge die ursprüngliche Matrix M
    M = create_fast_species_matrix(reactions, fast_species, slow_symbolic_variables, natnum, slow_species)

    # Modifiziere die Spalte S_prime in der Matrix M
    for S_double_prime in fast_species:
        sum_term = sum(
            multiply_slow_species(reaction, slow_symbolic_variables) *
            compute_scaled_rate(reaction, natnum, slow_species) *
            compute_species_difference(reaction, S)
            for reaction in reactions if consumes_fast_species(reaction, S_double_prime)
        )
        M[S_double_prime][S_prime] = sum_term  # Setze den neuen Wert für M[S_double_prime][S_prime]

    return M

def get_matrix_determinant(matrix):
    """
    Prüft, ob die gegebene Matrix eine Determinante von 0 hat bzw. gibt deren Determinante aus.
    
    :param matrix: Eine Dictionary-basierte symbolische Matrix {S: {T: Wert}}.
    :return: True, wenn die Determinante 0 ist, sonst False.
    """
    # Extrahiere die schnelle Spezies als sortierte Liste für eine feste Reihenfolge
    fast_species = sorted(matrix.keys())

    # Erstelle eine SymPy-Matrix aus der Dictionary-Struktur
    M_list = [[matrix[S][T] for T in fast_species] for S in fast_species]
    M_sympy = sp.Matrix(M_list)  # Konvertiere in SymPy-Matrix

    # Berechne die Determinante und vereinfache sie
    determinant = M_sympy.det().simplify()

    return determinant  # Prüft, ob die Determinante genau 0 ist

def sum_over_slow_reactions(data, natnum, slow_species, slow_symbolic_variables, slow_symbolic_derivatives):
    """
    Berechnet die erste Teilsumme des approximierten Generators für ein gegebenes Reaktionsnetzwerk.

    Diese Funktion summiert über alle langsamen Reaktionen, also Reaktionen, die nur langsame Spezies konsumieren 
    und produzieren. Für jede langsame Spezies wird ein Term zur Summe addiert, der sich aus folgenden Faktoren 
    zusammensetzt:
    - Der Differenz der Spezies in der Reaktion (compute_species_difference)
    - Der skalierten Reaktionsrate (compute_scaled_rate)
    - Dem Produkt der langsamen Spezies in der Reaktion (multiply_slow_species)
    - Der symbolischen Ableitung der langsamen Spezies (slow_symbolic_derivatives)

    :param data: Dictionary mit den Reaktionsdaten.
    :param natnum: Symbolische Variable für die Skalierung (z. B. N).
    :param slow_species: Liste der langsamen Spezies.
    :param slow_symbolic_variables: Symbolische Variablen der langsamen Spezies.
    :param slow_symbolic_derivatives: Symbolische Ableitungen der langsamen Spezies.
    :return: Die berechnete Summe als symbolischer Ausdruck.
    """
    total_sum = 0
    # Iteriere über jede Reaktion
    for reaction in data["reactions"]:
        # Prüfe, ob die Reaktion nur langsame Spezies konsumiert und produziert
        if is_slow_only_reaction(reaction, fast_species):
            # Summiere über alle langsamen Spezies
            for species in slow_species:
                species_diff = compute_species_difference(reaction, species)
                scaled_rate = compute_scaled_rate(reaction, natnum, slow_species)
                slow_species_product = multiply_slow_species(reaction, slow_symbolic_variables)
                species_derivative = slow_symbolic_derivatives[species]
                
                # Berechne das Produkt und addiere es zur Gesamtsumme
                total_sum += species_diff * scaled_rate * slow_species_product * species_derivative

    return total_sum

def sum_over_slow_fast_species_reactions(data, natnum, slow_species, fast_species, slow_symbolic_variables, slow_symbolic_derivatives):
    """
    Berechnet die zweite Teilsumme des approximierten Generators für das gegebene Reaktionsnetzwerk.

    Diese Funktion summiert über alle langsamen Spezies und schnellen Spezies und berücksichtigt dabei
    die Wechselwirkungen zwischen ihnen. Die zentrale Größe im Summanden ist der Ausdruck:

        det(M{slow, fast}) / det(M),

    wobei:
    - M die von `create_fast_species_matrix` erzeugte Matrix ist.
    - M{slow, fast} die von `create_modified_fast_species_matrix` erzeugte modifizierte Matrix ist, 
      in der die Spalte `fast` durch einen Vektor aus der b-Matrix ersetzt wurde.
    - Die Determinanten werden mit `get_matrix_determinant` berechnet.

    :param data: Dictionary mit den Reaktionsdaten.
    :param natnum: Symbolische Variable für die Skalierung (z. B. N).
    :param slow_species: Liste der langsamen Spezies.
    :param fast_species: Liste der schnellen Spezies.
    :param slow_symbolic_variables: Symbolische Variablen der langsamen Spezies.
    :param slow_symbolic_derivatives: Symbolische Ableitungen der langsamen Spezies.
    :return: Die berechnete Summe als symbolischer Ausdruck.
    """
    total_sum = 0

    # Berechne die ursprüngliche Matrix M einmal
    M = create_fast_species_matrix(data["reactions"], fast_species, slow_symbolic_variables, natnum, slow_species)
    det_M = get_matrix_determinant(M)
    
    # Iteriere über alle langsamen Spezies S
    for slow in slow_species:
        # Iteriere über alle schnellen Spezies S'
        for fast in fast_species:
            
            # Berechne die modifizierte Matrix M{slow, fast}
            M_modified = create_modified_fast_species_matrix(
                data["reactions"], fast_species, slow_symbolic_variables, natnum, slow_species, slow, fast
            )
            det_M_modified = get_matrix_determinant(M_modified)

            # Berechne den Ausdruck det(M{slow, fast}) / det(M)
            determinant_ratio = det_M_modified / det_M if det_M != 0 else 0 
            
            # Iteriere über alle Reaktionen
            for reaction in data["reactions"]:
                
                # Prüfe, ob die Reaktion die schnelle Spezies erzeugt
                if produces_fast_species(reaction, fast, fast_species):
                    
                    # Berechne die verschiedenen Terme
                    species_diff = compute_species_difference(reaction, slow)
                    scaled_rate = compute_scaled_rate(reaction, natnum, slow_species)
                    slow_species_product = multiply_slow_species(reaction, slow_symbolic_variables)
                    species_derivative = slow_symbolic_derivatives[slow]
                    
                    # Berechne den Summanden mit dem Determinantenverhältnis
                    summand = scaled_rate * slow_species_product * (species_diff + determinant_ratio) * species_derivative
                    
                    # Addiere zum Gesamtwert
                    total_sum += summand
    
    return total_sum

def total_sum_of_reactions(data, natnum, slow_species, fast_species, slow_symbolic_variables, slow_symbolic_derivatives):
    """
    Berechnet die Gesamtreaktionssumme für den approximierten Generator H^Nf.

    Diese Funktion berechnet die Gesamtreaktionssumme durch das Addieren zweier Teilsummen:
    1. Die Summe der langsamen Reaktionen (`slow_reactions_sum`).
    2. Die Summe der langsamen und schnellen Spezies-Reaktionen (`slow_fast_species_reactions_sum`).

    Nachdem beide Summen berechnet wurden, werden sie zusammenaddiert und der gesamte Ausdruck
    wird vereinfacht, um die endgültige Gesamtreaktionssumme zu erhalten.

    :param data: Dictionary mit den Reaktionsdaten.
    :param natnum: Symbolische Variable für die Skalierung (z. B. N).
    :param slow_species: Liste der langsamen Spezies.
    :param fast_species: Liste der schnellen Spezies.
    :param slow_symbolic_variables: Symbolische Variablen der langsamen Spezies.
    :param slow_symbolic_derivatives: Symbolische Ableitungen der langsamen Spezies.
    :return: Die vereinfacht berechnete Gesamtreaktionssumme.
    """
    # Berechne die Summe der langsamen Reaktionen
    slow_reactions_sum = sum_over_slow_reactions(data, natnum, slow_species, slow_symbolic_variables, slow_symbolic_derivatives)
    
    # Berechne die Summe der langsamen und schnellen Spezies-Reaktionen
    slow_fast_species_reactions_sum = sum_over_slow_fast_species_reactions(data, natnum, slow_species, fast_species, slow_symbolic_variables, slow_symbolic_derivatives)
    
    # Addiere beide Summen, um den approximierten Generator H^Nf zu berechnen
    total_sum = slow_reactions_sum + slow_fast_species_reactions_sum
    
    # Vereinfache den gesamten Ausdruck
    simplified_total_sum = simplify(total_sum)
    
    return simplified_total_sum








###################################




def is_crn_connected(reactions):
    # Graph als Adjazenzliste aufbauen
    graph = defaultdict(set)
    
    for reaction in reactions:
        educts = frozenset(reaction.get("educts", {}).keys())
        products = frozenset(reaction.get("products", {}).keys())
        
        if educts and products:
            graph[educts].add(products)
            graph[products].add(educts)
    
    # Falls der Graph leer ist, ist er nicht zusammenhängend
    if not graph:
        return False
    
    # Startknoten wählen (beliebiger erster Knoten im Graph)
    start = next(iter(graph))
    
    # DFS zur Überprüfung der Erreichbarkeit aller Knoten
    visited = set()
    stack = [start]
    
    while stack:
        node = stack.pop()
        if node not in visited:
            visited.add(node)
            stack.extend(graph[node] - visited)  # Nur unbesuchte Nachbarn hinzufügen
    
    # Der Graph ist zusammenhängend, wenn alle Knoten besucht wurden
    return len(visited) == len(graph)





def all_reactions_contain_fast_species(reactions, fast_species):
    fast_species_set = set(fast_species)  # Stelle sicher, dass fast_species ein Set ist
    for reaction in reactions:
        educts = set(reaction.get("educts", {}).keys())
        products = set(reaction.get("products", {}).keys())
        
        # Jetzt verwenden wir den Set-Operator für educts und products
        if not (educts & fast_species_set) or not (products & fast_species_set):
            return False
    return True





def find_connected_components(reactions):
    # Erstelle den Graphen als Adjazenzliste
    graph = defaultdict(set)
    reaction_mapping = {}
    
    for reaction in reactions:
        educts = frozenset(reaction.get("educts", {}).keys())
        products = frozenset(reaction.get("products", {}).keys())
        
        if educts and products:
            graph[educts].add(products)
            graph[products].add(educts)
            reaction_mapping[(educts, products)] = reaction
    
    # Falls der Graph leer ist, gibt es keine Komponenten
    if not graph:
        return []
    
    # Bestimme die verbundenen Komponenten
    visited = set()
    components = []
    
    for node in graph:
        if node not in visited:
            # Tiefensuche (DFS), um alle erreichbaren Knoten zu markieren
            stack = [node]
            connected_nodes = set()
            connected_reactions = []
            
            while stack:
                current = stack.pop()
                if current not in visited:
                    visited.add(current)
                    connected_nodes.add(current)
                    stack.extend(graph[current] - visited)
            
            # Finde zugehörige Reaktionen
            for (educts, products), reaction in reaction_mapping.items():
                if educts in connected_nodes or products in connected_nodes:
                    connected_reactions.append(reaction)
            
            components.append(connected_reactions)
    
    return components


def check_under_crns_for_fast_species(reactions, fast_species):
    # Schritt 1: Finde alle verbundenen Komponenten (Unter-CRNs)
    components = find_connected_components(reactions)
    
    # Schritt 2: Überprüfe für jede Unter-CRN, ob alle Reaktionen eine "fast" Spezies enthalten
    results = []
    for idx, component in enumerate(components):
        # Überprüfe die Unter-CRN (Komponente)
        if all_reactions_contain_fast_species(component, fast_species):
            results.append((idx, True))  # Unter-CRN idx enthält fast Spezies in allen Reaktionen
        else:
            results.append((idx, False))  # Unter-CRN idx enthält nicht in allen Reaktionen schnelle Spezies
    
    return results




###################################





def count_connected_components(reactions):
    # Graph als Adjazenzliste aufbauen
    graph = defaultdict(set)
    
    for reaction in reactions:
        educts = frozenset(reaction.get("educts", {}).keys())
        products = frozenset(reaction.get("products", {}).keys())
        
        if educts and products:
            graph[educts].add(products)
            graph[products].add(educts)
    
    # Falls der Graph leer ist, gibt es keine Komponenten
    if not graph:
        return 0
    
    # Anzahl der verbundenen Komponenten bestimmen
    visited = set()
    components = 0
    
    for node in graph:
        if node not in visited:
            # DFS zur Markierung aller erreichbaren Knoten
            stack = [node]
            while stack:
                current = stack.pop()
                if current not in visited:
                    visited.add(current)
                    stack.extend(graph[current] - visited)
            components += 1
    
    return components












###################################

def save_yaml(data, filename):
    """Speichert ein Dictionary als YAML-Datei."""
    with open(filename, "w") as file:
        yaml.dump(data, file, default_flow_style=False)

def extract_sub_crns(data):
    """Findet verbundene Komponenten und speichert sie als separate Sub-CRNs."""
    sub_crns = []
    components = find_connected_components(data["reactions"])  # Funktion existiert bereits
    
    for i, component in enumerate(components):
        sub_crn = {
            "name": f"{data['name']}_sub_{i+1}",
            "natnum": data["natnum"],
            "species": {"slow": [], "fast": []},
            "reactions": component
        }
        
        # Spezies extrahieren
        species_set = set()
        for reaction in component:
            species_set.update(reaction.get("educts", {}).keys())
            species_set.update(reaction.get("products", {}).keys())
        
        # Langsame und schnelle Spezies zuordnen
        sub_crn["species"]["slow"] = [s for s in data["species"]["slow"] if s in species_set]
        sub_crn["species"]["fast"] = [s for s in data["species"]["fast"] if s in species_set]
        
        sub_crns.append(sub_crn)
    
    return sub_crns

def save_sub_crns_as_yaml(data, filename):
    """Speichert die extrahierten Sub-CRNs als separate YAML-Dateien."""
    sub_crns = extract_sub_crns(data)
    base_name = filename.replace(".yaml", "")
    
    for i, sub_crn in enumerate(sub_crns):
        sub_filename = f"{base_name}sub_{i+1}.yaml"
        save_yaml(sub_crn, sub_filename)
        print(f"Saved: {sub_filename}")

###################################

def check_all_sub_crns_for_fast_species(reactions, fast_species):
    """
    Überprüft, ob in allen Sub-CRNs jede Reaktion mindestens eine schnelle Spezies
    sowohl in den Edukten als auch in den Produkten enthält.
    
    Args:
        reactions (list): Liste der Reaktionen im gesamten CRN.
        fast_species (list): Liste der schnellen Spezies.
    
    Returns:
        dict: Ein Dictionary mit den Sub-CRN-IDs als Schlüssel und `True` oder `False` als Wert.
    """
    # Schritt 1: Finde alle Sub-CRNs
    sub_crns = find_connected_components(reactions)

    # Schritt 2: Überprüfe für jedes Sub-CRN, ob alle Reaktionen die Bedingung erfüllen
    results = {}
    for idx, sub_crn in enumerate(sub_crns):
        results[f"sub_crn_{idx+1}"] = all_reactions_contain_fast_species(sub_crn, fast_species)

    return results





















#####################################



# Laden der Daten

#data = load_yaml("crn.yaml")

# Nutzer nach CRN-Datei fragen
filename = input("\nWhich CRN would you like to load? ") + ".yaml"
data = load_yaml(filename)

# Name des geladenen CRNs ausgeben
print(f"\nLoaded CRN: {data.get('name', 'Unbekannt')}")


slow_species = data["species"]["slow"]
fast_species = data["species"]["fast"]
natnum = sp.Symbol(data["natnum"])
slow_symbolic_variables = init_slow_symbolic_variables(slow_species)
slow_symbolic_derivatives = init_slow_symbolic_derivatives(slow_species)

#print(fast_species)
#print(slow_species)


print("\nIs this CRN connected? ", is_crn_connected(data["reactions"]))
print(f"\nNumber of connected sub-CRNs: {count_connected_components(data['reactions'])}")

#M = create_fast_species_matrix(data["reactions"], fast_species, slow_symbolic_variables, natnum, slow_species)
#det_M = get_matrix_determinant(M)

#print(det_M.simplify())


#results = check_under_crns_for_fast_species(data["reactions"], fast_species)

#print(results)






save_sub_crns_as_yaml(data, filename)


results = check_all_sub_crns_for_fast_species(data["reactions"], fast_species)

print(results)



toto = total_sum_of_reactions(data, natnum, slow_species, fast_species, slow_symbolic_variables, slow_symbolic_derivatives)

print("\nApproximate generator Hᴺf is given by")
sp.pprint(sympify(toto).simplify())  # Schöne symbolische Ausgabe




