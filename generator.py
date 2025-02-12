import yaml
import sympy as sp

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

data = load_yaml("crn.yaml")

slow_species = data["species"]["slow"]
natnum = sp.Symbol(data["natnum"])  # Natürliche Zahl als Symbol
slow_symbolic_variables = init_slow_symbolic_variables(slow_species)
slow_symbolic_derivatives = init_slow_symbolic_derivatives(slow_species)

reaction = data["reactions"][0]  # Beispiel: erste Reaktion
slow_educts = get_slow_educts(reaction, slow_species)

symbolic_product = multiply_slow_species(reaction, slow_symbolic_variables)
summed_coefficients = sum_slow_educt_coefficients(reaction, slow_species)
scaled_rate = compute_scaled_rate(reaction, natnum, slow_species)
species_difference = compute_species_difference(reaction, "S")

print(slow_educts)  # Erwartete Ausgabe: {'S': 44}
print(slow_symbolic_variables)  # Erwartete Ausgabe: {'S': v_{S}, 'PP': v_{PP}}
print(slow_symbolic_derivatives)  # Erwartete Ausgabe: {'S': df_{S}, 'PP': df_{PP}}
print(symbolic_product)  # Erwartete symbolische Ausgabe: v_{S}^{44}
print(summed_coefficients)  # Erwartete Ausgabe: 44
print(scaled_rate)  # Erwartete Ausgabe: k1 * N^(g1 - 1 + 44)
print(species_difference)  # Erwartete Ausgabe: Differenz aus Produkt- und Edukt-Koeffizienten für "S"
