from JSON_fill import load_json



def print_CRN(data):
    reactionstring = ""
    j = 0
    k = 0

    # Öffne eine Datei zum Schreiben (wenn die Datei existiert, wird sie überschrieben)
    with open("reactions_data_print.txt", "w") as file:
        
        while j < int(data["reactions_count"]):
            while k < int(data["slow_species_count"]):
                if int(data["reactions"][j]["educts"][f"S{k+1}"]) == 0:
                    k = k + 1
                else:
                    reactionstring = reactionstring + "+ " + str(data["reactions"][j]["educts"][f"S{k+1}"]) + " * S" + str(k+1) + " "
                    k = k + 1
            k = 0
            
            Ej_index = None  # Variable für den gefundenen Index

            for key, value in data["reactions"][j]["educts"].items():
                if key.startswith("E") and value == 1:  # Suche nach "E" mit Wert 1
                    Ej_index = key[1:]  # Extrahiere die Zahl hinter "E"
                    break  # Es gibt genau eine "E" mit Wert 1, also abbrechen

            # Erzeuge den reactionstring basierend auf dem gefundenen Index
            if Ej_index is None:
                reactionstring += " ---> "
            else:
                reactionstring += f"+ E{Ej_index} ---> "

            while k < int(data["slow_species_count"]):
                if data["reactions"][j]["products"][f"S{k+1}"] == 0:
                    k = k + 1
                else:
                    reactionstring = reactionstring + "+ " + str(int(data["reactions"][j]["products"][f"S{k+1}"])) + " * S" + str(k+1) + " "
                    k = k + 1
            k = 0

            Ej_index = None  # Variable für den gefundenen Index

            for key, value in data["reactions"][j]["products"].items():
                if key.startswith("E") and value == 1:  # Suche nach "E" mit Wert 1
                    Ej_index = key[1:]  # Extrahiere die Zahl hinter "E"
                    break  # Es gibt genau eine "E" mit Wert 1, also abbrechen

            # Erzeuge den reactionstring basierend auf dem gefundenen Index
            if Ej_index is None:
                reactionstring += ""
            else:
                reactionstring += f"+ E{Ej_index}"
            
            # Schreibe den reactionstring und die anderen Daten in die Datei
            file.write(f"Reaktion {j+1}:\n")
            file.write(reactionstring + "\n")
            file.write("κ_ℛ = " + data["reactions"][j]["rate_unscaled"] + ", γ_ℛ = " + data["reactions"][j]["scaling"] + "\n"+ "\n")
            reactionstring = str()
            j = j + 1
        j = 0
        file.write("\n")  # Leere Zeile am Ende der Datei

    
    
# Lade die JSON-Daten von der Datei reaction_data.json
data = load_json("reaction_data.json")

# Ausgabe des CRN
print_CRN(data)