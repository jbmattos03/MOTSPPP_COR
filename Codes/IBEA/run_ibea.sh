#!/usr/bin/env bash
set -euo pipefail

echo "Compilando o programa IBEA MOTSPPP..."
g++ -O3 -o ibea_motsppp index.cpp structures.cpp BoundedParetoSet.cpp || { echo "Erro na compilação."; exit 1; }

sizes=(5 8 10 20 50 100 200)

for t in "${sizes[@]}"; do
  # instâncias simétricas
  for i in {1..6}; do
    dir="../resultados/ibea/a3/sym/${t}.${i}"
    echo "Executando instância ${t}.${i} simétrica..."
    rm -rf "$dir"
    mkdir -p "$dir"
    ./ibea_motsppp "../instances/A3/symmetric/${t}.${i}.in" "${dir}/dummy.txt"
  done

  # instâncias assimétricas
  for i in {1..6}; do
    dir="../resultados/ibea/a3/asym/${t}.${i}"
    echo "Executando instância ${t}.${i} assimétrica..."
    rm -rf "$dir"
    mkdir -p "$dir"
    ./ibea_motsppp "../instances/A3/asymmetric/${t}.${i}.in" "${dir}/dummy.txt"
  done
done

rm -f ibea_motsppp
echo "Todas as execuções foram concluídas!"