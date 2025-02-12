import json
import sys


def load_json(filename="reaction_data.json"):
    with open(filename, "r") as f:
        return json.load(f)

def validate_json(data):
    errors = []  # Liste zum Sammeln von Fehlern

    # Schritt 1: Schleife über Reaktionen
    for reaction in data["reactions"]:
        reaction_index = reaction["index"]

        # Schritt 2: Überprüfen der educts S Keys
        for key, value in reaction["educts"].items():
            if key.startswith("S"):
                if not (value == "" or (isinstance(value, str) and value.isdigit() and int(value) >= 0) or (isinstance(value, int) and value >= 0)):
                    errors.append(f"❌ Invalid value in educts[{key}] for reaction {reaction_index}: '{value}'. Allowed values: '', 'x' (where x >= 0), or natural number.")
        
        # Schritt 3: Überprüfen der educts E Keys
        e_count = 0  # Zähler für "1" oder 1 in E Keys
        for key, value in reaction["educts"].items():
            if key.startswith("E"):
                if value not in ["", "0", "1", 0, 1]:  # Werte, die erlaubt sind
                    errors.append(f"❌ Invalid value in educts[{key}] for reaction {reaction_index}: '{value}'. Allowed values: '', '0', '1'.")
                if value in ["1", 1]:
                    e_count += 1
                    if e_count > 1:
                        errors.append(f"❌ More than one '1' or 1 found in educts[{key}] for reaction {reaction_index}. Only one '1' or 1 is allowed.")

        # Schritt 4: Überprüfen der products S Keys
        for key, value in reaction["products"].items():
            if key.startswith("S"):
                if not (value == "" or (isinstance(value, str) and value.isdigit() and int(value) >= 0) or (isinstance(value, int) and value >= 0)):
                    errors.append(f"❌ Invalid value in products[{key}] for reaction {reaction_index}: '{value}'. Allowed values: '', 'x' (where x >= 0), or natural number.")

        # Schritt 5: Überprüfen der products E Keys
        e_count = 0  # Zähler für "1" oder 1 in E Keys
        for key, value in reaction["products"].items():
            if key.startswith("E"):
                if value not in ["", "0", "1", 0, 1]:  # Werte, die erlaubt sind
                    errors.append(f"❌ Invalid value in products[{key}] for reaction {reaction_index}: '{value}'. Allowed values: '', '0', '1'.")
                if value in ["1", 1]:
                    e_count += 1
                    if e_count > 1:
                        errors.append(f"❌ More than one '1' or 1 found in products[{key}] for reaction {reaction_index}. Only one '1' or 1 is allowed.")

        # Schritt 6: Überprüfen von rate_unscaled
        #if not isinstance(reaction["rate_unscaled"], str):
            #errors.append(f"❌ Invalid value in rate_unscaled for reaction {reaction_index}. Should have been a string. Run this program a second time.")

        # Schritt 7: Überprüfen von scaling
        #if not isinstance(reaction["scaling"], str):
            #errors.append(f"❌ Invalid value in scaling for reaction {reaction_index}. Should have been a string. Run this program a second time.")

    if errors:
        # Alle Fehler anzeigen
        for error in errors:
            print(error)
        sys.exit("❌ Program terminated due to validation errors.")
    else:
        print("✅ All checks passed successfully!")

def fix_values(data):
    # Schritt 8: Fixieren der Werte ("" auf 0 und "x" auf 0 für natürliche Zahlen)
    for reaction in data["reactions"]:
    # Umwandeln in Edicts
        for key, value in reaction["educts"].items():
            if key.startswith("S") or key.startswith("E"):
                if value == "":
                    reaction["educts"][key] = 0
                elif isinstance(value, str) and value.isdigit():
                    reaction["educts"][key] = int(value)
                elif isinstance(value, float):  # Falls ein Float vorhanden ist, auf Integer setzen
                    reaction["educts"][key] = int(value)

    # Umwandeln in Products
        for key, value in reaction["products"].items():
            if key.startswith("S") or key.startswith("E"):
                if value == "":
                    reaction["products"][key] = 0
                elif isinstance(value, str) and value.isdigit():
                    reaction["products"][key] = int(value)
                elif isinstance(value, float):  # Falls ein Float vorhanden ist, auf Integer setzen
                    reaction["products"][key] = int(value)
                    
    
    # rate_unscaled und scaling anpassen
        if isinstance(reaction["rate_unscaled"], int):  # Umwandeln in String
            reaction["rate_unscaled"] = str(reaction["rate_unscaled"])
        if reaction["rate_unscaled"] == "":
                reaction["rate_unscaled"] = f"k{reaction['index']}"

        if isinstance(reaction["scaling"], int):  # Umwandeln in String
            reaction["scaling"] = str(reaction["scaling"])
        if reaction["scaling"] == "":
                reaction["scaling"] = f"g{reaction['index']}"

    print("✅ Fixed values for 'S' and 'E' keys: Empty strings set to 0, strings with numbers converted to integers.")
    print("✅ Fixed values for 'rate_unscaled' and 'scaling' keys: Empty strings set to reaction-index-variables, numbers converted to strings.")

def save_to_json(data, filename="reaction_data.json"):
    with open(filename, "w") as f:
        json.dump(data, f, indent=4)
    print(f"Data successfully saved to {filename}")

# Lade die JSON-Daten von der Datei reaction_data.json
data = load_json("reaction_data.json")

# Führe die Fehlerüberprüfung durch
validate_json(data)

# Korrigiere die Werte
fix_values(data)

# Speichere die geänderten Daten zurück in die Datei
save_to_json(data)
