c     programa liasp2
c     calcula a dissimilaridade entre uma rede e aquelas obtidas pela elimina��o
c     de cada um de seus n�s individualmente. ignora contribui��o do n� eliminado
c     matrizes de vizinhanca calculada com madchar13
c     entra dados no formato de uma unica matriz de adjacencia
c     onde cada linha � uma sequencia sem espa�os de 0's e 1's
c     ex: 0111111110111100010101000000000100

      parameter(npm=1000)
      integer*1 am(npm,npm),bm(npm,npm),lar(2,npm)
      integer*2 mv1(npm*(npm-1)/2),mv2(npm*(npm-1)/2),kki(0:npm,0:npm)
      integer*2 mv3(npm*(npm-1)/2)
      character*100 entrada3,saida5

      open (unit=2,file='liasp2.dat')

c     Parametros e variaveis auxiliares:
c     npm = numero maximo de nos
c     am = matriz de adjacencia (utilizado na manipulacao da rede)
c     bm = matriz de adjacencia (copia para nao perder a rede original)
c     lar = vetor auxiliar para leitura da matriz de adjacencia
c     mv1 = vetor de vizinhanca
c     mv2 = vetor de vizinhanca
c     mv3 = vetor de vizinhanca
c     kki = ??
c     entrada3 = nome do arquivo de entrada
c     saida5 = nome do arquivo de saida

c     Variaveis de entrada:
c     nm = numero maximo de nos na rede
c     np = potencia booleana maxima

c     Variaveis de saida:
c     nc = numero de nos conectados na rede
c     na = numero de arestas na rede
c     mcl = tamanho do maior cluster
c     nmcl = numero de maiores cluster
c     diam = diametro do maior cluster
c     xmed = dissimilaridade m�dia entre os n�s do maior cluster
c     dissimilaridade = dissimilaridade entre a rede original e aquela com o no eliminado


c     entrada de dados
c------------------------------------------------------------------
c o arquivo de entrada tem a seguinte estrutura:
c     nm np [num maximo de nós na rede, potencia booleana maxima]
c     entrada3 [nome do arquivo de entrada com matriz de adjacencia]
c     saida5 [nome do arquivo de saida]
c==================================================================

c     leitura dos parametros
10    read (2,*)nm,np
      if(nm.lt.0)stop
      read (2,500)entrada3
      read (2,500)saida5

c     abre arquivo da matriz de adjacencia
      open (unit=3,file=entrada3)

c     prepara o arquivo de saida
      open (unit=5,file=saida5)
      write(5,*)'#  no elim  nc     na     mcl    nmcl   diam    xmed
     .    dist1       dist2'

c     inicializa o vetor auxiliar para leitura da matriz de adjacencia
      do i = 1,npm
       lar(1,i) = 0
      enddo

c------------------------------------------------------------------
c     le a matriz de adjacencia
c------------------------------------------------------------------

c     inicializa contadores de numero de arestas [na] e nos conectados [nc]
      na = 0
      nc = 0

c     le a matriz de adjacencia, linha a linha
c     armazenando-a em um vetor auxiliar [lar] e contando o numero de nos conectados
c     e o numero de arestas
c     ic = indicador de no conectado (1 se o no esta conectado, 0 caso contrario)
      do i = 1,nm
       ic = 0
       read (3,510)(lar(1,j),j=1,nm)
       do j = 1,nm
        am(i,j) = lar(1,j)
        bm(i,j) = lar(1,j)
        na = na + am(i,j)
        ic = max(ic,am(i,j))
       enddo
       nc = nc + ic
      enddo

c     chama o programa madchar13 para calcular a matriz de vizinhanca
c     da rede original, armazenando-a em mv1
      call madchar13(am,mv1,nm,nc,np,kki,xmed,id,mcl,idege)

c     fecha o arquivo da matriz de adjacencia
      ii = 0
      d = 0.
      write(5,520)ii,nc,na,mcl,idege,id,xmed,d,d

c==================================================================
c     comeco do grand loop de eliminacao de nos
c==================================================================

c     para cada no da rede...
      do ii = 1,nm

c     reverte a rede para a rede original
       do i = 1,nm
        do j = 1,nm
         am(i,j) = bm(i,j)
        enddo
       enddo

c     elimina o no ii da rede, zerando a linha e a coluna ii
       do j = 1,nm
        am(ii,j) = 0
        am(j,ii) = 0
       enddo

c     calcula o numero de nos conectados [nc] e o numero de arestas [na]
c     (semelhante ao que foi feito na leitura da matriz de adjacencia)
       na = 0
       nc = 0
       do i = 1,nm
        ic = 0
        do j = 1,nm
         na = na + am(i,j)
         ic = max(ic,am(i,j))
        enddo
        nc = nc + ic
       enddo

c     chama o programa madchar13 para calcular a matriz de vizinhanca da rede
c     com o no ii eliminado e armazena o resultado em mv2
       call madchar13(am,mv2,nm,nc,np,kki,xmed,id,mcl,idege)

c==================================================================
c     gera mv3 zerando os elementos de mv1 correspondentes ao no eliminado
c     isto e util para separar a contribuicao direta e indireta do no eliminado
c     para o calculo da dissimilaridade.

c     faz uma copia de mv1 em mv3
       do i = 1,nm*(nm-1)/2
        mv3(i) = mv1(i)
       enddo

c     inicializa o acumulador de dissimilaridade direta
       dissim_d = 0.

c     zera os elementos de mv3 correspondentes ao no ii e acumula a dissimilaridade
c     direta
       do i = 1,nm-1
        do j = i+1,nm
         if(i.eq.ii.or.j.eq.ii)then
          ll = (2*nm-i)*(i-1)/2+j-i
c         acumula a dissimilaridade direta
c         a dissimilaridade direta e a diferenca entre a eficiencia
c         da rede original e a eficiencia da rede com o no ii eliminado
c         que ocorre exclusivamente devido as arestas que ligam o no ii
c         aos seus vizinhos
          if(mv3(ll).ne.0)then
           dissim_d = dissim_d + 2./mv3(ll)
          mv3(ll) = 0
         endif
        enddo
       enddo

c     divide a dissimilaridade direta acumulada pelo numero de arestas
c     possiveis na rede original
       dissim_d = dissim_d/(nm*(nm-1))

c==================================================================
c     calcula a dissimilaridade total entre mv1 e mv2
c     (dissimilaridade entre as matrizes de vizinhanca)
       dissim_t = dissim(mv1,mv2,nm)

c      escreve os resultados no arquivo de saida
       write(5,520)ii,nc,na,mcl,idege,id,xmed,dissim_t,dissim_d

c     fim do grand loop de eliminacao de nos
      enddo

c     fecha o arquivo de saida
      close (unit=5)

c     retorna ao inicio do programa
c     (o programa termina quando chega ao fim do arquivo de entrada,
c     indicado pela presenca de um 'numero de nos' negativo)
      goto 10

c     declaracao dos formatos de saida utilizados
500   format(a100)
510   format(10000i1)
520   format(6(1x,i6),3(2x,e10.3))
540   format(10000i4)

      stop

      end

c     fim do programa principal
c==================================================================
c==================================================================

      subroutine madchar13(a,mv,nm,nc,np,kk,xlmd,id,mcl,idege)
c==================================================================

c     Parametros e variaveis auxiliares:
c     lla = ??
c     pro = ??
c     pra = ??
c     lis = lista auxiliar de conexoes
c           armazena qual no esta conectado a qual de forma compacta.
c           Uma matriz de tamanho (numero max de arestas na rede) x 2.
c           O i-esimo elemento da segunda linha da matriz indica o indice
c           do elemento da primeira linha no qual comeca a lista de nos
c           conectados ao no i.
c           Os elementos da primeira linha armazenados entre os indices lis(i,2)
c           e lis(i+1,2) sao os nos conectados ao no i.
c     ga = ??
c     kmd = ??
c     gamd = ??
c     lmd = ??
c     set = ??
c     nmcl = numero de maiores cluster
c     diam = diametro do maior cluster
c     xmed = dissimilaridade media entre os nos do maior cluster
c     dissimilaridade = dissimilaridade entre a rede original e aquela com o no eliminado

c     Argumentos
c     1: a = matriz de adjacencia (entrada)
c     2: mv = matriz de vizinhanca (saida)
c     3: nm = numero de nos na rede (entrada)
c     4: nc = numero de nos conectados (entrada)
c     5: np = potencia booleana maxima (entrada)
c     6: kk = ?? (saida)
c     7: xlmd = dissimilaridade media entre os nos do maior cluster (saida)
c     8: id = diametro da rede (saida)
c     9: mcl = tamanho do maior cluster (saida)
c     10: idege = numero de arestas da rede (saida)


      parameter(npm=1000)
      integer*1 a(npm,npm)
      integer*2 mv(npm*(npm-1)/2),kk(0:npm,0:npm),lla(npm),pro(npm)
      integer*2 pra(npm)
      integer lis((npm-1)*npm/2,2)
      real ga(0:npm), kmd(0:npm)
      real gamd(0:npm),lmd(0:npm),set(npm,npm)
c==================================================================
c     coloca zero no valores medios das dissimilaridades entre nos
      xlmd = 0.
      ylmd = 0.
      zlmd = 0.

      do i = 1,npm
       lla(i) = 0
       do j = 1,npm
        set(i,j) = 0.
       enddo
      enddo

      do i = 1,npm*(npm-1)/2
       mv(i) = 0
      enddo

      do i = 0,npm
       lmd(i) = 0.
       ga(i) = 0
       kmd(i) = 0
       gamd(i) = 0
       do j = 0,npm
        kk(i,j) = 0
       enddo
      enddo
c==================================================================
c     escreve a(i,j) em entrada
c     varre a triagular superior da matriz de adjacencia e preenche
c     os valores da triangular inferior.
c     armazena os valores da triangular superior em um vetor mv.
c     a variavel auxiliar ll converte os indices da triangular superior
c     em indices do vetor mv.

      do i = 1,nm-1
       do j = i+1,nm
        a(j,i) = a(i,j)
        ll = (2*nm-i)*(i-1)/2+j-i
        mv(ll) = a(i,j)
       enddo
      enddo

c    define a diagonal principal da matriz de adjacencia como 1
c    pra? pro?
      do i = 1,nm
        a(i,i) = 1
        pra(i) = 1
        pro(i) = 0
      enddo
c==================================================================
c     prepara a lista de elementos nao nulos de a(i,j)

c     inicializa lis com zeros
      do i = 1,nm
       do j = 1,2
        lis(i,j) = 0.
       enddo
      enddo
c     iq e um contador de quantas arestas ja foram armazenadas em lis.
c     ele e utilizado para indicar a posicao em que foi armazenada a
c     lista de nos conectados a a cada no.
      iq = 1

c     para cada linha da matriz de adjacencia (cada no da rede)
c     armazena o indice da linha em que comeca a lista de nos conectados
c     ao no i
      do i = 1,nm
       lis(i,2) = iq
       do j = 1,nm
        if(a(i,j).eq.1)then
c        para cada no conectado ao no i, armazena o indice do no
c        na lista de nos conectados ao no i e incrementa o contador iq
         lis(iq,1) = j
         iq = iq + 1
        endif
       enddo
      enddo

c    ao sair do loop, a variavel i vale nm+1.
c    o elemento lis(nm+1,2) indica o indice do ultimo elemento de lis(:, 1).
      lis(i,2) = iq
      
c==================================================================
c     comeca o loop para calculo as propriedades das diferentes matrizes mad(ip)

      do i = 1,nm-1
       do ip = 1,np
        if(pra(i).eq.0.and.pro(i).lt.ip)goto 90

        if (ip.eq.np) go to 100
c==================================================================
c     comeca a obtencao de am**ip - le mad e transfere para am
c==================================================================
c     le a matriz tmad5 faz pb tmad5*mad e guarda resultado em tmad6
c==================================================================

        pra(i) = 0

        do j = i+1,nm

         a(i,j) = 0

         if (a(j,i).eq.1) goto 126

         do k = lis(j,2),lis(j+1,2)-1
          k1 = lis(k,1)
          if(i.lt.k1) then
           if (mv((2*nm-i)*(i-1)/2+k1-i).gt.0)then
            if(mv((2*nm-i)*(i-1)/2+k1-i).le.ip)go to 124
           endif
          else if (i.gt.k1) then
           if (mv((2*nm-k1)*(k1-1)/2+i-k1).gt.0)then
            if(mv((2*nm-k1)*(k1-1)/2+i-k1).le.ip)go to 124
           endif
          endif
         enddo

         goto 125
124      a(i,j) = 1
         pra(i) = 1
125      continue
         ll = (2*nm-i)*(i-1)/2+j-i
         mv(ll) = a(i,j)*(ip+1)
         pro(j) = max(pro(j),mv(ll))

126     enddo

c==================================================================
c     obtem a nova matriz com vizinhos exclusivos de ordem ip+1: am=tmad6-tmad5
c==================================================================
c     transfere a nova matriz com todos os vizinhos de tdma1 para tdma0
c==================================================================
c     atualiza a matriz de adjacencia ponderada, le tmad1, soma (ip+1)*am e coloca em tmad2
c==================================================================
c     transfere a nova matriz com todos os vizinhos ponderados de tmad2 para tmad1
c==================================================================

        do j = i+1,nm
         a(j,i) = a(j,i) + a(i,j)
        enddo

       enddo

c======================================================
c     fim do loop em ip

90    enddo
c======================================================
c     saida dos resultados para a amostra is

100   npp = 1

      do i = 1,nm-1
       do j = i+1,nm
        ll = (2*nm-i)*(i-1)/2+j-i
        npp = max(npp,mv(ll))
       enddo
      enddo

c======================================================
c     calcula a ordem de sub-grafo de cada no

      do i = 1,nm
c       a(i,i) = 0
       do j = i+1,nm
        a(i,j) = a(j,i)
       enddo
      enddo


      do i = 1,nm
       kk(i,0) = 0
       do j = 1,nm
        kk(i,0) = kk(i,0) + a(i,j)
       enddo
      enddo

      do i = 1,nm
c       a(i,i) = 0
       kk(i,i) = 0
      enddo

      ll = 0
      do i = 1,nm-1
       do j = i+1,nm
        ll = ll + 1
c        a(i,j) = mv(ll)
c        a(j,i) = mv(ll)
        kk(i,j) = mv(ll)
        kk(j,i) = mv(ll)
       enddo
      enddo

      do i = 1,nm
       lmd(i) = 0
       do j = 1,nm
c        lmd(i) = lmd(i) + a(i,j)
        lmd(i) = lmd(i) + kk(i,j)
       enddo
      enddo
c======================================================
c     calcula a numero de nos no maior cluster

      mcl = 0
      do i = 1,nm
       mcl = max(kk(i,0),mcl)
      enddo
c======================================================
c     calcula a dissimilaridade minima m�dia da rede limitada ao(s) maior(es) cluster(s)

c     nm2 = (0*nm+nc)*(0*nm+nc)
      xm2 = 0.
      if (mcl.gt.1)xm2 = 1./float((mcl-1)*mcl)
      imc = 0
      do i = 1,nm
       if (kk(i,0).eq.mcl.and.kk(i,0)*nc.gt.0)then
        imc = imc + 1
        xlmd = xlmd + lmd(i)
c       xlmd = xlmd + lmd(i)*xm2
        ylmd = ylmd + lmd(i)/kk(i,0)/kk(i,0)
        zlmd = zlmd + lmd(i)/kk(i,0)/(0*nm+nc)
       endif
      enddo
      idege = imc/mcl
c      write(5,*)xlmd,xm2,mcl,imc
       xlmd = xlmd*xm2*mcl/max(imc,1)
c      write(5,*)xlmd,xm2,mcl,imc
c======================================================
c     transfere a matriz ponderada para o campo mv(np,np)
c======================================================
c     calcula diametro

      id = 0
      do i = 1,nm*(nm-1)/2
       id = max(mv(i),id)
      enddo
c======================================================
c     calcula dimensao fractal
c======================================================
c     escreve dissimilaridade minima media de cada no em saida9
c======================================================
c     escreve coeficiente de clusteriza��o em saida7
c======================================================
c     escreve coeficientes de assortatividade em saida8
c======================================================
c     escreve distribui��o de nos em saida8
c=====================================================
      return

      end

c     fim da subroutine madchar13(a,mv,nm,np,xmd,id)
c=====================================================================
c=====================================================================

      subroutine rede1(a,n,np,kk,ga,ip)
c=====================================================================
      parameter(npm=1000)
      integer*1 a(npm,npm)
      integer*2 kk(0:npm,0:npm)
      real ga(0:npm)
      integer vz(np),ordem
c=====================================================================
c     calcula o grau de cada no

      kk(0,ip) = 0
      do i = 1,n
       kk(i,ip) = 0
       do j = 1,n
        if(i.ne.j)then
         ii = min(i,j)
         jj = max(i,j)
         kk(i,ip) = kk(i,ip) + a(ii,jj)
        endif
       enddo
       kk(0,ip) = kk(0,ip) + kk(i,ip)
      enddo

      if (ip.gt.0) return

c=====================================================================
c     calcula o coeficiente de aglomera��o de cada noh

      if (ip.gt.1)return

      ga(0) = 0.
      do i=1,n
       do j=1,n
        vz(j) = -1
       enddo

       ordem = 0
       do j=1,n
        if (a(i,j) .eq. 1) then
         ordem=ordem+1
         vz(ordem) = j
        endif
       enddo

       ctotal = ((ordem * (ordem-1))/2)
       sumg=0
       do j=1,ordem
        li = vz(j)
        do k=1,n
         if (a(li,k).eq.1.and.a(i,k).eq.1.and.li.ne.k.and.k.ne.i)then
          sumg=sumg+1
         endif
        enddo
       enddo

       clocal = (sumg)/2
       if (ctotal .gt. 0) then
        ga(i) = clocal/ctotal
       else
        ga(i)=0
       endif
       ga(0) = ga(0) + ga(i)
      enddo

      return
      end

c     fim da subroutine rede1(a,n,np,nvm,kk,ga,ip)
c=====================================================================
      real function dissim(mv1,mv2,nm)
c     calcula a dissimilaridade entre duas matrizes de vizinhanca
      parameter(npm=1000)
      integer*2 mv1(npm*(npm-1)/2),mv2(npm*(npm-1)/2)

      dissim = 0.

      do i = 1,nm*(nm-1)/2

c       Se mv2(i) e' 0, o aumento no acumulador de dissimilaridade vale
c       1/mv1(i). Se ambos sao 0, nao ha' aumento. Caso ambos sejam
c       diferentes de 0, o aumento vale (1/mv1(i) - 1/mv2(i)).
       if(mv2(i) == 0)then
        if(mv1(i) > 0)then
         dissim = dissim +  2*1/mv1(i)
        endif
       elseif(mv1(i) > 0)then
c       So' checando se mv2(i) >= mv1(i):
        if(mv2(i) < mv1(i))then
         print(5,*)"! Erro: mv2(i) < mv1(i) !"
        endif
c ------------------------------------------
        dissim = dissim + 2*(1/mv1(i) - 1/mv2(i))
       endif
      enddo

      dissim = dissim/(nm-1)/nm

      return
      end

c     fim da real function dissim(mv1,mv2,nm)
c=====================================================================
c=====================================================================
