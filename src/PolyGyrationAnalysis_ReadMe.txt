# Como usar esse programa ?
A - Essa função foi modificada para realizar os cálculos de todos os frames ou vc pode passar um númerero específico de frames
Nesse exemplo estamos selecionando 1000 frames

    r = calc_rg("box_20.pdb", "9.1_prod_reduced.dcd",num_frames=1000, seg_select="segname POLI")

Nesse segundo exemplo não estamos passando nenhum número de frames.

    r = calc_rg("box_20.pdb", "9.1_prod_reduced.dcd", seg_select="segname POLI")

Os componentes dessa função:

    r = calc_rg("Arquivo_pdb - (obrigatória)", "trajetória - "obrigatória"", num_frames - (opcional)  seg_select="selecão do polimero - (obrigatória)")

B - Essa segunda etapa o plot do dado vai ser gerado

    plot_rg(r, "radius_10") 

Nesse caso a função vai receber dois argumentos. O primeiro vai ser os dados gerados na etapa anterior e o nome do arquivo de saida.

C - Por fim o scrip pode salvar a analise com a função:

    sv_rg(r, "radius_10.dat")


COMANDOS RÁPIDOS
    r = calc_rg("box_20.pdb", "9.1_prod_reduced.dcd",num_frames=1000, seg_select="segname POLI")
    plot_rg(r, "radius_10")
    sv_rg(r, "radius_10.dat")

