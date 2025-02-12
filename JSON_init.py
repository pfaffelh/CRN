import json

def get_input():
    data = {}

    # Eingabe für reactions_count (nur echte positive ganze Zahl)
    while True:
        try:
            data["reactions_count"] = int(input("How many reactions?    "))
            if data["reactions_count"] <= 0:
                print("❌ Invalid input. Please enter a positive integer (1, 2, 3, ...).")
            else:
                break
        except ValueError:
            print("❌ Invalid input. Please enter a positive integer (1, 2, 3, ...).")

    # Eingabe für slow_species_count (natürliche Zahl, 0 erlaubt)
    while True:
        try:
            data["slow_species_count"] = int(input("How many slow species? "))
            if data["slow_species_count"] < 0:
                print("❌ Invalid input. Please enter a non-negative integer (0, 1, 2, ...).")
            else:
                break
        except ValueError:
            print("❌ Invalid input. Please enter a non-negative integer (0, 1, 2, ...).")

    # Eingabe für fast_species_count (natürliche Zahl, 0 erlaubt)
    while True:
        try:
            data["fast_species_count"] = int(input("How many fast species? "))
            if data["fast_species_count"] < 0:
                print("❌ Invalid input. Please enter a non-negative integer (0, 1, 2, ...).")
            else:
                break
        except ValueError:
            print("❌ Invalid input. Please enter a non-negative integer (0, 1, 2, ...).")

    reactions = []

    for j in range(data["reactions_count"]):
        reaction = {}
        reaction["index"] = j + 1
        reaction["educts"] = {}  # Verwende ein Dictionary
        reaction["products"] = {}  # Verwende ein Dictionary

        # Slow species als educts
        for k in range(data["slow_species_count"]):
            reaction["educts"][f"S{k+1}"] = ""  # Leeres Feld, keine Zuweisung

        # Fast species als educts
        for k in range(data["fast_species_count"]):
            reaction["educts"][f"E{k+1}"] = ""  # Leeres Feld, keine Zuweisung

        # Slow species als products
        for k in range(data["slow_species_count"]):
            reaction["products"][f"S{k+1}"] = ""  # Leeres Feld, keine Zuweisung

        # Fast species als products
        for k in range(data["fast_species_count"]):
            reaction["products"][f"E{k+1}"] = ""  # Leeres Feld, keine Zuweisung

        reaction["rate_unscaled"] = ""  
        reaction["scaling"] = ""  
        
        reactions.append(reaction)

    data["reactions"] = reactions

    return data

def save_to_json(data, filename="reaction_data.json"):
    with open(filename, "w") as f:
        json.dump(data, f, indent=4)
    print(f"Data successfully saved to {filename}")

if __name__ == "__main__":
    reaction_data = get_input()
    save_to_json(reaction_data)
