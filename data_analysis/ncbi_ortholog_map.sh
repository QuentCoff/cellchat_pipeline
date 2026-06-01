#!/bin/bash

API_KEY="b5a7f4376319caa754def6debbda68da1f09"
INPUT="/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/data_analysis/gene_ids_list.txt"
OUTPUT="/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/data_analysis/ortholog_mapping_test.csv"

# Créer le fichier CSV avec en-tête s'il n'existe pas
if [ ! -f "$OUTPUT" ]; then
  echo "horse_gene_id,human_gene_id,human_symbol,status" > "$OUTPUT"
fi

# Liste des gènes déjà traités
already_done=$(cut -d',' -f1 "$OUTPUT" | tail -n +2)

total=$(wc -l < "$INPUT")
i=0
skipped=0

while read gid; do
  i=$((i+1))

  # Skip si déjà traité
  if echo "$already_done" | grep -qx "$gid"; then
    skipped=$((skipped + 1))
    continue
  fi

  echo ""
  echo ">>> [$i/$total] Gène cheval ID: $gid (skippés: $skipped)"

  # Retry loop (max 3 tentatives)
  attempt=0
  http_code="000"
  body=""
  while [ $attempt -lt 3 ]; do
    attempt=$((attempt + 1))
    resp=$(curl -s -w "\n%{http_code}" \
      --max-time 30 \
      -X GET "https://api.ncbi.nlm.nih.gov/datasets/v2/gene/id/${gid}/orthologs?returned_content=COMPLETE&taxon_filter=9606&page_size=10" \
      -H "accept: application/json" \
      -H "api-key: ${API_KEY}")

    http_code=$(echo "$resp" | tail -1)
    body=$(echo "$resp" | sed '$d')

    if [ "$http_code" = "200" ]; then
      break
    fi

    echo "    Tentative $attempt/3 échouée (HTTP $http_code), retry dans 2s..."
    sleep 2
  done

  if [ "$http_code" != "200" ]; then
    echo "    ERREUR HTTP $http_code après 3 tentatives"
    echo "$gid,,,HTTP_${http_code}" >> "$OUTPUT"
    sleep 1
    continue
  fi

  if [ -z "$body" ]; then
    echo "    ERREUR: réponse vide"
    echo "$gid,,,EMPTY_BODY" >> "$OUTPUT"
    sleep 0.2
    continue
  fi

  result=$(echo "$body" | python3 -c "
import sys, json
try:
    data = json.load(sys.stdin)
    reports = data.get('reports', [])
    if not reports:
        print('NO_ORTHOLOG')
    else:
        r = reports[0]
        g = r.get('gene', {})
        tax = g.get('tax_id', '')
        if str(tax) == '9606':
            hid = g.get('gene_id', 'NA')
            sym = g.get('symbol', 'NA')
            print(f'{hid},{sym}')
        else:
            print('NO_ORTHOLOG')
except Exception as e:
    print(f'PARSE_ERROR:{e}')
")

  if echo "$result" | grep -q "NO_ORTHOLOG"; then
    echo "    Pas d'orthologue humain"
    echo "$gid,,,NO_ORTHOLOG" >> "$OUTPUT"
  elif echo "$result" | grep -q "PARSE_ERROR"; then
    echo "    ERREUR parsing: $result"
    echo "$gid,,,${result}" >> "$OUTPUT"
  else
    echo "    → Humain: $result"
    echo "$gid,${result},OK" >> "$OUTPUT"
  fi

  # Délai anti-rate-limit (3 req/s max)
  sleep 0.35

done < "$INPUT"

echo ""
echo ">>> TERMINE. Résultat dans: $OUTPUT"
