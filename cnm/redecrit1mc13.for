c     programa redecrit1mc13 
c     igual a redecrit1mc11, mas usa a versao madchar13 para o calculo da matriz de vizinhanca 
c     calcula distancia entre duas redes determinadas para valores distintos 
c     de uma vari�vel de controle 
c     a partir de suas matrizes de vizinhanca calculada com madchar13 
c     entra dados no formato de uma matriz de intera��o dependente do periodo 
c     unica diferen�a � escrever a matriz de vizinhan�a para o �ltimo valor de ip
c     entra com variaveis reais, ao passo que redecrit0mc11 entra variaveis inteiras
c     vers�o para tratar redes do projeto epigen

      parameter(npm=1500)
      integer*1 am(npm,npm)
      integer*2 mv1(npm*(npm-1)/2),mv2(npm*(npm-1)/2),kki(0:npm,0:npm)
      real lar(2,npm)
      real idis(npm,npm)
      character*100 entrada3,saida5,saida6,saida7

      open (unit=2,file='redecrit1mc13.dat')

c     Parametros e variaveis auxiliares:
c     npm = numero maximo de nos
c     am = matriz de adjacencia (utilizado na manipulacao da rede)
c     mv1 = matriz de vizinhanca em formato de vetor (triangular superior)
c     mv2 = matriz de vizinhanca em formato de vetor (triangular superior)
c     mv3 = matriz de vizinhanca em formato de vetor (triangular superior)
c     kki = ??
c     idis = matriz de distancias
c     xp = treshold do parametro de controle
c     entrada3 = nome do arquivo de entrada
c     saida5 = nome do arquivo de saida (propriedades da rede para cada valor de xp)
c     saida6 = nome do arquivo de saida (ultima matriz de vizinhanca)
c     saida7 = nome do arquivo de saida (ultima matriz de adjacencia)

c     Variaveis de entrada:
c     nm = numero maximo de nos na rede
c     np = potencia booleana maxima
c     pmin = valor minimo do parametro
c     ns = numero de valores do parametro
c     delp = intervalo de variacao do parametro

c     Variaveis de saida:
c     nc = numero de nos conectados na rede
c     na = numero de arestas na rede
c     mcl = numero de nos no maior cluster
c     xmed = dissimilaridade media entre os nos da rede filtrada

c     Colunas do arquivo de saida5

c     1: xp = valor do parametro de controle para o qual a rede foi filtrada
c     2: nc = numero de nos conectados na rede
c     3: na = numero de arestas na rede original
c     4: mcl = numero de nos no maior cluster
c     5: idege = numero de arestas na rede filtrada
c     6: id = diametro da rede filtrada
c     7: xmed = dissimilaridade media entre os nos da rede filtrada
c     8: d(ii+1) = dissimilaridade entre as redes ii e ii+1

c     entrada de dados
c------------------------------------------------------------------
c     o arquivo de entrada tem a seguinte estrutura:
c      nm np [num maximo de nós na rede, potencia booleana maxima]
c      pmin, delp, ns [valor minimo, intervalo de variacao do parametro, numero de valores de xp]
c      entrada3 [nome do arquivo de entrada com matriz de adjacencia]
c      saida5 [nome do arquivo de saida]
c      saida6 [nome do arquivo de saida]
c      saida7 [nome do arquivo de saida]
c==================================================================

c     leitura dos parametros
10    read (2,*)nm,np
      if(nm.lt.0)stop
      read (2,*)pmin,delp,ns
      read (2,500)entrada3
      read (2,500)saida5
      read (2,500)saida6
      read (2,500)saida7

c     abre arquivo da matriz de conexoes (e.g. correlation matrix)
      open (unit=3,file=entrada3)

c     prepara os arquivos de saida
      open (unit=5,file=saida5)
      open (unit=6,file=saida6)
      open (unit=7,file=saida7)
      write(5,*)'#      rede   nc   na  mcl nmcl diam        xmed    d(i,
     .i+1)'
      write(6,*)'#      rede   nc   na  mcl nmcl diam        xmed    d(i,
     .i+1)'

c     inicializa o vetor auxiliar para leitura da matriz de conexoes
      do i = 1,npm
       lar(1,i) = 0
      enddo

c------------------------------------------------------------------
c     le a matriz de conexoes
c------------------------------------------------------------------
c     le a matriz de conexoes, linha a linha
c     armazenando-a em um vetor auxiliar [lar]

      do i = 1,nm
       read (3,*)(lar(1,j),j=1,nm)
       do j = 1,nm
        idis(i,j) = lar(1,j)
       enddo
c      define a diagonal principal como zero
       idis(i,i) = 0.
      enddo

c     inicializa o treshold do parametro de controle [xp]
c     esse valor comeca com o valor minimo [pmin] e vai sendo incrementado
c     de acordo com o intervalo [delp]
      xp = pmin

c     Transforma a matriz de conexoes em uma matriz de adjacencia
c     a matriz de adjacencia [am] e uma matriz binaria, onde
c     am(i,j) = 1 se a conexao entre os nos i e j for maior que o treshold,
c     am(i,j) = 0 caso contrario.
c     inicializa contadores de numero de arestas [na] e nos conectados [nc]
c     ic = indicador de no conectado (1 se o no esta conectado, 0 caso contrario)
      na = 0
      nc = 0
      do i = 1,nm
       ic = 0
       do j = 1,nm
        am(i,j) = 0
c         if(idis(i,j).le.xp.and.idis(i,j).gt.0)then 
         if(idis(i,j).ge.xp.and.idis(i,j).gt.0)then 
         am(i,j) = 1
         na = na + 1
         ic = 1
        endif
       enddo
       nc = nc + ic
      enddo
c     divide o numero de arestas por 2, pois cada aresta e contada duas vezes
      na = na/2

c     chama o programa madchar13 para calcular a matriz de vizinhanca
      call madchar13(am,mv1,nm,nc,np,kki,xmed,id,mcl,idege)

c     escreve as propriedades iniciais da rede no arquivo de saida 5
c     inicializa a variavel de varredura do loop de amostras [is]
c     a distancia entre da rede original aa ela mesma e zero
      is = 1
      d = 0.
      write(5,520)xp,nc,na,mcl,idege,id,xmed,d

c================================================================== 
c     comeco do grand loop no numero de amostras
c==================================================================

      do is = 1,ns

c     atualiza o valor do treshold
       xp = pmin + is*delp

c     Transforma a matriz de conexoes em uma matriz de adjacencia
c     utilizando o treshold atualizado.
c     (similar ao que foi feito acima)
       na = 0
       nc = 0
       do i = 1,nm
        ic = 0
        do j = 1,nm
         am(i,j) = 0
c         if(idis(i,j).le.xp.and.idis(i,j).gt.0)then 
         if(idis(i,j).ge.xp.and.idis(i,j).gt.0)then
          am(i,j) = 1
          na = na + 1
          ic = 1
         endif
        enddo
        nc = nc + ic
       enddo
       na = na/2

c     caso o numero de arestas nao seja zero...
       if (na.gt.0) then

c     chama o programa madchar13 para calcular a matriz de vizinhanca
c     atualizada [mv2] e calcula a distancia entre as duas redes
c     [mv1] e [mv2]
        call madchar13(am,mv2,nm,nc,np,kki,xmed,id,mcl,idege)
        d = dist(mv1,mv2,nm)

       else

c     caso o numero de arestas seja zero, a distancia entre as duas redes
c     e definida como 0.
        na = 0
        mcl = 0
        xmed = 0
        idege = 0
        id = 0
        d = 0.

       endif

c     escreve os resultados no arquivo de saida 5
       write(5,520)xp,nc,na,mcl,idege,id,xmed,d

c     atualiza a matriz de vizinhanca [mv1] para a proxima iteracao
       do i = 1,nm*(nm-1)/2
         mv1(i) = mv2(i)
       enddo
      enddo

c     Transforma a matriz de adjacencia do formato vetorial (elementos da
c     diagonal superior) para o formato matricial (matriz simetrica)
c     utilizando a variavel auxiliar idis.
      ll = 0
      do i = 1,nm
       do j = i+1,nm
        ll = ll + 1
        idis(i,j) = mv1(ll)
        idis(j,i) = mv1(ll)
       enddo
      enddo

c    escreve a matriz de vizinhanca no arquivo de saida 6
      do i = 1,nm
       write(6,540)(int(idis(i,j)),j=1,nm)
      enddo

c     atualiza a matriz de adjacencia com base na matriz de vizinhanca.
c     Se a distancia entre dois nos nao for 1, eles nao sao adjacentes.
      do i = 1,nm
       do j = 1,nm
        if(idis(i,j).ne.1)am(i,j)=0
       enddo
      enddo

c     escreve a matriz de adjacencia no arquivo de saida 7
      do i = 1,nm
       write(7,510)(am(i,j),j=1,nm)
      enddo

c     fecha os arquivos de saida
      close (unit=3)
      close (unit=5)
      close (unit=6)
      close (unit=7)

c     retorna ao inicio do programa
c     (o programa termina quando chega ao fim do arquivo de entrada,
c     indicado pela presenca de um 'numero de nos' negativo)
      goto 10

c     declaracao dos formatos de saida utilizados
500   format(a100) 
510   format(10000i1) 
520   format(1x,f10.5,5(1x,i6),2(2x,e10.3)) 
c520  format(6(1x,i4),2(2x,e10.3))
c530   format(10000(f5.0,1x)) 
c530   format(10000(i6,1x)) 
540   format(10000i4)

      stop

      end

c     fim do programa principal
c==================================================================
c==================================================================

      subroutine madchar13(a,mv,nm,nc,np,kk,xlmd,id,mcl,idege)
c==================================================================

c     Parametros e variaveis auxiliares:
c     pro = maior caminho entre dois nos da rede
c     pra = vetor de indicadores.
c           inicializado como 1, ao entrar no loop de calculo de caminhos
c           a partir do no' i, pra(i) e' alterado para 0. Caso haja um
c           caminho entre o no' i e o no' j, pra(j) e' alterado para 0.
c     lis = lista auxiliar de conexoes
c           armazena qual no esta conectado a qual de forma compacta.
c           Uma matriz de tamanho (numero max de arestas na rede) x 2.
c           O i-esimo elemento da segunda linha da matriz indica o indice
c           do elemento da primeira linha no qual comeca a lista de nos
c           conectados ao no i.
c           Os elementos da primeira linha armazenados entre os indices lis(i,2)
c           e lis(i+1,2) sao os nos conectados ao no i.
c     lmd = vetor com a soma dos caminhos minimos que ligam o no' i a 
c           todos os outros nos.

c     Argumentos
c     1: a = matriz de adjacencia (entrada)
c     2: mv = matriz de vizinhanca (saida)
c     3: nm = numero de nos na rede (entrada)
c     4: nc = numero de nos conectados (entrada)
c     5: np = potencia booleana maxima (entrada)
c     6: kk = matriz com a ordem de subgrafo de cada no armazenada na coluna 0,
c             os demais elementos da matriz armazenam a matriz de vizinhanca. (saida)
c     7: xlmd = comprimento medio dos menores caminhos entre nos da rede
c               limitada ao(s) maior(es) subgrafo(s) conexo(s) (saida)
c     8: id = diametro da rede (saida)
c     9: mcl = tamanho do maior cluster (saida)
c     10: idege = numero de arestas da rede (saida)


      parameter(npm=1000)
      integer*1 a(npm,npm)
      integer*2 mv(npm*(npm-1)/2),kk(0:npm,0:npm), pro(npm)
      integer*2 pra(npm)
      integer lis((npm-1)*npm/2,2)
      real ga(0:npm), kmd(0:npm)
      real lmd(0:npm)
c==================================================================
c     inicializa as variaveis
      xlmd = 0.
c      ylmd = 0. (soma dos caminhos minimos normalizada pelo quadrado
c                 do numero de nos no maior subgrafo conexo)
c      zlmd = 0. (soma dos caminhos minimos normalizados pelo produto
c                 do numero de nos no maior subgrafo conexo e o numero
c                 de nos na rede)

      do i = 1,npm*(npm-1)/2
       mv(i) = 0
      enddo

      do i = 0,npm
       lmd(i) = 0.
       kmd(i) = 0
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

c     para cada no da rede...
c     ( loop em i )
      do i = 1,nm-1 
      
c     ip e' o contador de potencia booleana (?)
c     ( loop em ip )
       do ip = 1,np

c       ?????
        if(pra(i).eq.0.and.pro(i).lt.ip)goto 90

c       caso ip chegue ao maior caminho a ser considerado, saia do loop
        if (ip.eq.np) go to 100

c       atualiza pra(i), indicando que estao sendo calculados os caminhos
c       que ligam o no' i a todos os outros nos da rede
        pra(i) = 0

c       para cada no, varre os indices da triangular superior da matriz de adjacencia
c       ( loop em j )
        do j = i+1,nm

c        zera o valor da conexao entre os nos i e j
         a(i,j) = 0

c        caso os nos i e j estejam conectados, 
c        a informacao ja esta armazenada.
c        passa para o proximo no'.
         if (a(j,i).eq.1) goto 126

c        varre os nos conectados ao no j, armazenados em lis(:,1).
c        os vizinhos de j estao armazenados entre os indices lis(j,2) e lis(j+1,2)-1.
c        ( loop em k )
         do k = lis(j,2), lis(j+1,2)-1

c         armazena o indice do no vizinho do no j em k1
          k1 = lis(k,1)

c         verifica se existe um caminho que liga o no' k1 (vizinho de j) ao no' i com tamanho
c         menor ou igual a ip. se sim, entao os nos i e j estao conectados.
c         Armazena um caminho com tamanho ip+1 entre i e j no vetor da 
c         matriz de vizinhanca [mv]. A condicional se i < k1 e' necessaria
c         pois mv so' armazena a triangular superior.
          if(i.lt.k1) then
           if (mv((2*nm-i)*(i-1)/2+k1-i).gt.0)then
            if(mv((2*nm-i)*(i-1)/2+k1-i).le.ip)go to 124
           endif
          else if (i.gt.k1) then
           if (mv((2*nm-k1)*(k1-1)/2+i-k1).gt.0)then
            if(mv((2*nm-k1)*(k1-1)/2+i-k1).le.ip)go to 124
           endif
          endif

c        ( final do loop em k )
         enddo

c        pula a atribuicao de 1 a a(i,j) e vai para o proximo no j
         goto 125

c        caso exista caminho ligando os nos i e j:
124      a(i,j) = 1
         pra(i) = 1

c        armazena o menor caminho na matriz de vizinhanca
125      continue
         ll = (2*nm-i)*(i-1)/2+j-i
         mv(ll) = a(i,j)*(ip+1)
c        pro(j) armazena o maior caminho minimo ja encontrado com uma
c        extremidade no no j.
         pro(j) = max(pro(j), mv(ll))

c       ( final do loop em j )
126     enddo

c       atualiza a triangular inferior da matriz a, indicando quais
c       nos estao conectados ao no i, com um caminho de tamanho arbitrario.
c       a triangular inferior armazena informacao se o existe caminho ja 
c       encontrado entre os nos i e j, com tamanho arbitrario.
c       a triangular superior armazena informacao se foi encontrado um caminho
c       entre os nos i e j para o tamanho de caminho ip+1, que esta sendo calculado.
        do j = i+1,nm
         a(j,i) = a(j,i) + a(i,j)
        enddo

c       ( final do loop em ip )
       enddo

c     ( final do loop em i )
90    enddo
100   continue

c     npp nao e' utilizado
c100   npp = 1
c
c      do i = 1,nm-1
c       do j = i+1,nm
c        ll = (2*nm-i)*(i-1)/2+j-i
c        npp = max(npp, mv(ll))
c       enddo
c      enddo

c======================================================
c     calcula a ordem de sub-grafo de cada no'

c     define a triangular superior da matriz de adjacencia
c     igual aa triangular inferior. A matriz a armazena informacao
c     se existe um caminho entre os nos i e j, com tamanho arbitrario.
c     Dessa forma, e' possivel identificar os subgrafos de cada no.
      do i = 1,nm
       do j = i+1,nm
        a(i,j) = a(j,i)
       enddo
      enddo

c     atualiza a coluna 0 da matriz kk, que armazena a ordem de sub-grafo
c     de cada no' correspondente a cada linha.
      do i = 1,nm
       kk(i,0) = 0
       do j = 1,nm
        kk(i,0) = kk(i,0) + a(i,j)
       enddo
      enddo

c     define a diagonal principal da matriz kk como 0
      do i = 1,nm
       kk(i,i) = 0
      enddo

c     atualiza os demais elementos da matriz kk, que armazena a matriz de vizinhanca
      ll = 0
      do i = 1,nm-1
       do j = i+1,nm
        ll = ll + 1
        kk(i,j) = mv(ll)
        kk(j,i) = mv(ll)
       enddo
      enddo

c     calcula a soma dos caminhos minimos que ligam o no' i a todos os outros nos
      do i = 1,nm
       lmd(i) = 0
       do j = 1,nm
        lmd(i) = lmd(i) + kk(i,j)
       enddo
      enddo

c     calcula a numero de nos no maior subgrafo conexo.
c     O numero de nos no subgrafo do qual o no' i faz parte
c     e' igual a kk(i,0).
      mcl = 0
      do i = 1,nm
       mcl = max(kk(i,0), mcl)
      enddo
c======================================================
c     calcula o comprimento medio dos menores caminhos entre nos da rede
c     limitada ao(s) maior(es) subgrafo(s) conexo(s)

c     inicializa e armazena o inverso do numero maximo de arestas
c     possiveis no maior subgrafo conexo.
      xm2 = 0.
      if (mcl.gt.1)xm2 = 1./float((mcl-1)*mcl)

c     inicializa um contador de numero de nos
      imc = 0

c     calcula 
      do i = 1,nm
       if (kk(i,0).eq.mcl.and.kk(i,0)*nc.gt.0)then
        imc = imc + 1
        xlmd = xlmd + lmd(i)
c       xlmd = xlmd + lmd(i)*xm2
c       ylmd = ylmd + lmd(i)/kk(i,0)/kk(i,0)
c       zlmd = zlmd + lmd(i)/kk(i,0)/(0*nm+nc)
       endif
      enddo

c     idege e' a razao entre o numero de nos considerado no calculo
c     da dissimilaridade e o numero total de nos.
c     pode ser entendido como o numero de subgrafos conexos com o
c     numero maximo de nos.
c     idege > 1. indica que ha' mais de um subgrafo conexo com o mesmo
c     numero de nos.
      idege = imc/mcl

c     calcula 
      xlmd = xlmd*xm2*mcl/max(imc,1)

c      write(5,*)xlmd,xm2,mcl,imc

c     calcula diametro da rede (diametro do maior subgrafo conexo)
      id = 0
      do i = 1,nm*(nm-1)/2
       id = max(mv(i),id)
      enddo

      return

      end

c     fim da subroutine madchar13(a,mv,nm,np,xmd,id)
c=====================================================================
c=====================================================================

      subroutine rede1(a,n,np,kk,ga,ip)
c=====================================================================
      parameter(npm=1500)
      integer*1 a(npm,npm)
      integer*2 kk(0:npm,0:npm)
      real ga(0:npm)
      integer vz(np)
      integer ordem
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
      real function dist(mv1,mv2,nm)
c     calcula a dissimilaridade entre duas matrizes de vizinhanca
      parameter(npm=1000)
      integer*2 mv1(npm*(npm-1)/2),mv2(npm*(npm-1)/2)

      dist = 0.

      do i = 1,nm*(nm-1)/2

c       Se mv2(i) e' 0, o aumento no acumulador de dissimilaridade vale
c       1/mv1(i). Se ambos sao 0, nao ha' aumento. Caso ambos sejam
c       diferentes de 0, o aumento vale (1/mv1(i) - 1/mv2(i)).
       if(mv2(i) == 0)then
        if(mv1(i) > 0)then
         dist = dist +  2*1/mv1(i)
        endif
       elseif(mv1(i) > 0)then
c       So' checando se mv2(i) >= mv1(i):
        if(mv2(i) < mv1(i))then
         print *, "! Erro: mv2(i) < mv1(i) !"
        endif
c ------------------------------------------
        dist = dist + 2*(1/mv1(i) - 1/mv2(i))
       endif
      enddo

      dist = dist/(nm-1)/nm

      return
      end


c     fim da real function dist(mv1,mv2,nm)
c=====================================================================
c=====================================================================
