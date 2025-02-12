import yaml
import json

# YAML-Datei laden
with open("network.yaml", "r") as f:
    yaml_data = yaml.safe_load(f)

# YAML-Daten als JSON speichern
with open("network.json", "w") as f:
    json.dump(yaml_data, f, indent=4)

print("✅ YAML wurde erfolgreich in JSON umgewandelt!")
