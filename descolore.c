#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

/*------------------------------------------------
 * Fila de prioridade
 *------------------------------------------------*/
typedef struct no *ptno;

typedef struct no {
    int i, j;
    ptno next;
} no;

/*-----------------------------------------------
 *  init   *Q      new     init         *Q
 *  [a:]->[b:]-+   [c:]  ->[a:]->[b:]->[c:]-+
 *   ^---------+             ^--------------+
 *----------------------------------------------*/
void insQ(ptno *Q, int i, int j) {
    ptno new = malloc(sizeof (no));
    new->i = i;
    new->j = j;
    if (!(*Q))
        new->next = new;
    else {
        new->next = (*Q)->next;
        (*Q)->next = new;
    }
    *Q = new;
}

/*-----------------------------------------------
 *  init         *Q                    *Q
 *  [a:]->[b:]->[c:]-+    [a:]  [b:]->[c:]-+
 *    ^--------------+            ^--------+
 *----------------------------------------------*/
void remQ(ptno *Q, int *i, int *j) {
    ptno init = (*Q)->next;
    *i = init->i;
    *j = init->j;
    if (*Q == init)
        *Q = NULL;
    else
        (*Q)->next = init->next;
    free(init);
}

int isEmpty(ptno *Q) {
    return *Q == NULL;
}

void initQPrior(ptno *QPrior, int mn) {
    int i;
    for (i = 0; i < mn; i++)
        QPrior[i] = NULL;
}

void insert(ptno *QPrior, int i, int j, int p) {
    insQ(QPrior + p, i, j);
}

int pop(ptno *QPrior, int *i, int *j, int *maxPrior, int mn) {
    while (*maxPrior < mn && isEmpty(QPrior + *maxPrior))
        (*maxPrior)++;
    if (*maxPrior == mn)
        return 1;
    remQ(QPrior + *maxPrior, i, j);
    return 0;
}

/*------------------------------------------------
 * Imagem PGM
 *------------------------------------------------*/
#define CREATOR "# Imagem criada por Luiz Eduardo da Silva\n"

typedef int *image;

image img_alloc(int nl, int nc) {
    return (image) malloc(nl * nc * sizeof (int));
}

int img_free(image Img) {
    free(Img);
}

void img_info(char *nome, int nl, int nc, int mn) {
    printf("\nInformacoes:");
    printf("\n--------------------------\n");
    printf("Nome do arquivo.............: %s \n", nome);
    printf("Numero de linhas............: %d \n", nl);
    printf("Numero de colunas...........: %d \n", nc);
    printf("Maximo nivel-de-cinza/cores.: %d \n\n", mn);
}

void ERROR(char *s) {
    printf("ERRO: %s\n", s);
    exit(10);
}

image read_pgm(char *nome, int *nl, int *nc, int *mn) {
    int i, j, k;
    char c, LINHA[100];
    image Img;
    FILE *arq;
    if ((arq = fopen(nome, "r")) == NULL)
        ERROR("Abertura do arquivo PGM");
    //- PGM = "P2" -----------/
    fgets(LINHA, 80, arq);
    if (!strstr(LINHA, "P2"))
        ERROR("Formato do arquivo PGM");
    //- le as linhas de comentario --/
    do {
        fgets(LINHA, 80, arq);
    } while (LINHA[0] == '#');
    //- le as dimensoes da imagem --/
    sscanf(LINHA, "%d %d", nc, nl);
    fscanf(arq, "%d", mn);

    if (*nl == 0 || *nc == 0 || *mn == 0)
        ERROR("Dimensões da imagem PGM");
    Img = img_alloc(*nl, *nc);
    if (!Img)
        ERROR("Alocação de memória para imagem");
    for (i = 0; i < *nl * (*nc); i++) {
        fscanf(arq, "%d", &k);
        if (k > *mn)
            ERROR("valores dos pixels na imagem");
        Img[i] = k;
    }
    fclose(arq);
    return Img;
}

void read_ppm(image *R, image *G, image *B, char *nome, int *nl, int *nc, int *mn) {
    int i, j, r, g, b;
    char LINHA[100];
    FILE *arq;
    if ((arq = fopen(nome, "r")) == NULL)
        ERROR("Abertura do arquivo PGM");
    /*-- PGM = "P2" -----------*/
    fgets(LINHA, 80, arq);
    if (!strstr(LINHA, "P3"))
        ERROR("Formato do arquivo PGM");
    /*-- le as linhas de comentario --*/
    do {
        fgets(LINHA, 80, arq);
    } while (LINHA[0] == '#');
    /*-- le as dimensoes da imagem --*/
    sscanf(LINHA, "%d %d", nc, nl);
    fscanf(arq, "%d", mn);

    if (*nl == 0 || *nc == 0 || *mn == 0)
        ERROR("Dimensões da imagem PGM");

    //printf("%d %d %d", nl, nc, mn);

    *R = img_alloc(*nl, *nc);
    *G = img_alloc(*nl, *nc);
    *B = img_alloc(*nl, *nc);

    if (*R && *G && *B) {
        for (i = 0; i < *nl; i++) {
            for (j = 0; j < *nc; j++) {
                fscanf(arq, "%d %d %d", &r, &g, &b);
                if (r > *mn || g > *mn || b > *mn)
                    ERROR("Erro nos dados do arquivo \n\nn");
                (*R)[i * (*nc) + j] = r;
                (*G)[i * (*nc) + j] = g;
                (*B)[i * (*nc) + j] = b;
            }
        }
        fclose(arq);
    }
    else
        ERROR("Erro na alocacao de memoria para o arquivo \n\n");
}

void write_pgm(image Img, char *nome, int nl, int nc, int mn) {
    int i, j, x, k, valores_por_linha;
    FILE *arq;
    if ((arq = fopen(nome, "wt")) == NULL)
        ERROR("Criação do arquivo PGM");
    fputs("P2\n", arq);
    fputs(CREATOR, arq);
    fprintf(arq, "%d  %d\n", nc, nl);
    fprintf(arq, "%d\n", mn);
    valores_por_linha = 20;
    for (i = 0, k = 0; i < nl; i++)
        for (j = 0; j < nc; j++) {
            x = Img[i * nc + j];
            fprintf(arq, "%3d ", x);
            k++;
            if (k > valores_por_linha) {
                fprintf(arq, "\n");
                k = 0;
            }
        }
    fclose(arq);
}

void write_ppm(image R, image G, image B, char *nome, int nl, int nc, int mn) {
    int i, j, r, g, b, k, valores_por_linha = 16;
    FILE *arq;
    if ((arq = fopen(nome, "wt")) == NULL)
        ERROR("Erro na criacao do arquivo");
    fputs("P3\n", arq);
    fputs(CREATOR, arq);
    fprintf(arq, "%d  %d\n", nc, nl);
    fprintf(arq, "%d\n", mn);
    k = 0;
    for (i = 0; i < nl; i++)
        for (j = 0; j < nc; j++) {
            r = R[i * nc + j];
            g = G[i * nc + j];
            b = B[i * nc + j];
            fprintf(arq, "%3d %3d %3d ", r, g, b);
            k++;
            if (k > valores_por_linha) {
                fprintf(arq, "\n");
                k = 0;
            }
        }
    fclose(arq);
}

double min(double R, double G, double B) {
    if (R <= G && R < B)
        return R;
    else if (G <= R && G < B)
        return G;
    else
        return B;
}

void rbg2hsi(double R, double G, double B, double *H, double *S, double *I) {
    double num, den, theta;

    num = .5 * ((R - G) + (R - B));
    den = pow((R - G) * (R - G) + (R - B) * (G - B), .5);
    theta = acos(num / (den + 0.00001)) / (2 * M_PI) * 360;

    *H = (B > G) ? 360 - theta : theta;
    *S = 1 - 3.0 / (R + G + B + 0.00001) * min(R, G, B);
    *I = 1 / 3.0 * (R + G + B + 0.00001);
}

void hsi2rgb(double H, double S, double I, double *R, double *G, double *B) {
    if (H < 120) {
        *B = I * (1 - S);
        *R = I * (1 + S * cos(H / 360 * 2 * M_PI) / cos((60 - H) / 360 * 2 * M_PI));
        *G = 3 * I - (*R + *B);
    }
    else if (H < 240) {
        H = H - 120;
        *R = I * (1 - S);
        *G = I * (1 + S * cos(H / 360 * 2 * M_PI) / cos((60 - H) / 360 * 2 * M_PI));
        *B = 3 * I - (*R + *G);
    }
    else {
        H = H - 240;
        *G = I * (1 - S);
        *B = I * (1 + S * cos(H / 360 * 2 * M_PI) / cos((60 - H) / 360 * 2 * M_PI));
        *R = 3 * I - (*G + *B);
    }
}

void gradient(image In, image Out, int nl, int nc, int mn, int raio) {
    int i, j, y, x, max, min;
    for (i = 0; i < nl; i++)
        for (j = 0; j < nc; j++) {
            max = -1;
            min = mn + 1;
            for (y = -raio; y <= raio; y++) // -1 0 1
                for (x = -raio; x <= raio; x++) // -1 0 1
                {
                    int pi = i + y;
                    int pj = j + x;
                    if (pi >= 0 && pi < nl && pj >= 0 && pj < nc) {
                        if (abs(x) + abs(y) <= raio && In[pi * nc + pj] > max)
                            max = In[pi * nc + pj];
                        if (abs(x) + abs(y) <= raio && In[pi * nc + pj] < min)
                            min = In[pi * nc + pj];
                    }
                }
            Out[i * nc + j] = max - min;
        }
}

void watershed(image In, image Out, int nl, int nc, int mn, int y, int x) {
    int i, j, k, maxPrior = 0, stop = 0;
    ptno qPrior[mn];
    image mark = img_alloc(nl, nc);

    enum {
        NONE,
        QUEUE,
        WSHED,
        MARK1,
        MARK2,
    };

    struct {
        int i, j;
    } n4[4] = {0, 1, 1, 0, 0, -1, -1, 0};

    initQPrior(qPrior, mn);
    // Inicialização dos marcadores
    //-----------------------------
    for (i = 0; i < nl * nc; i++)
        mark[i] = NONE;
    // MARK1 = (y, x) - centro da região

    int raio = 10;
    for (i = -raio; i <= raio; i++)
        for (j = -raio; j <= raio; j++) {
            int pi = i + y;
            int pj = j + x;
            if (abs(i) + abs(j) <= raio)
                mark[pi * nc + pj] = MARK1;
        }

    //mark[y * nc + x] = MARK1;

    // MARK2 = borda da imagem
    for (i = 0; i < nl; i++) {
        mark[i * nc] = MARK2;
        mark[i * nc + nc - 1] = MARK2;
    }
    for (j = 0; j < nc; j++) {
        mark[j] = MARK2;
        mark[(nl - 1) * nc + j] = MARK2;
    }
    // Inicialização
    //--------------
    // Todos os vizinhos dos marcadores são colocados na fila
    // A prioridade é o nível de cinza do pixel
    for (i = 1; i < nl - 1; i++)
        for (j = 1; j < nc - 1; j++)
            if (mark[i * nc + j] == NONE) {
                int isAdj = 0;
                for (k = 0; k < 4; k++) {
                    int pi = i + n4[k].i;
                    int pj = j + n4[k].j;
                    int m = mark[pi * nc + pj];
                    if (m == MARK1 || m == MARK2)
                        isAdj = 1;
                }
                if (isAdj) {
                    mark[i * nc + j] = QUEUE;
                    insert(&qPrior, i, j, In[i * nc + j]);
                }
            }
    // Crescimento dos Marcadores
    //---------------------------
    // Enquanto a fila de prioridade não está vazia faça
    // 1. retira um pixel da fila
    // 2. Se o pixel eh vizinho de UMA marca, recebe essa marca.
    //    Seus vizinhos sem rotulo são colocados na fila
    // 3. Se o pixel eh vizinho de DUAS marcas, ele eh WATERSHED
    while (!stop) {
        int m = NONE;
        int isWshed = 0;
        stop = pop(qPrior, &i, &j, &maxPrior, mn);
        if (!stop) {
            for (k = 0; k < 4; k++) {
                int pi = i + n4[k].i;
                int pj = j + n4[k].j;
                if (pi >= 0 && pi < nl && pj >= 0 && pj < nc) {
                    int mAdj = mark[pi * nc + pj];
                    if (mAdj == MARK1 || mAdj == MARK2) {
                        if (m == NONE)
                            m = mAdj;
                        else if (m != mAdj)
                            isWshed = 1;
                    }
                }
            }
            if (isWshed)
                mark[i * nc + j] = WSHED;
            else {
                mark[i * nc + j] = m;
                for (k = 0; k < 4; k++) {
                    int pi = i + n4[k].i;
                    int pj = j + n4[k].j;
                    if (pi >= 0 && pi < nl && pj >= 0 && pj < nc)
                        if (mark[pi * nc + pj] == NONE) {
                            int prior, px;
                            mark[pi * nc + pj] = QUEUE;
                            px = In[pi * nc + pj];
                            prior = (px < maxPrior) ? maxPrior : px;
                            insert(&qPrior, pi, pj, prior);
                        }
                }
            }
        }
    }
    for (i = 0; i < nl * nc; i++)
        Out[i] = (mark[i] == MARK2) ? 255 : 0;
    img_free(mark);
}

void tonsCinzaParc(image R, image G, image B, image Out, int nl, int nc, int mn) {
    int i;
    double H, S, I, a, b, c;

    for (i = 0; i < nl * nc; i++) {
        if (Out[i] == 0) {
            a = (double) R[i] / mn;
            b = (double) G[i] / mn;
            c = (double) B[i] / mn;

            rbg2hsi(a, b, c, &H, &S, &I);

            // Zerar H e S para deixar em tons de cinza
            H = 0;
            S = 0;

            //transforma HSI >> RGB
            hsi2rgb(H, S, I, &a, &b, &c);
            R[i] = a * mn;
            G[i] = b * mn;
            B[i] = c * mn;
        }
    }
}

void tonsCinza(image R, image G, image B, int nl, int nc, int mn) {
    int i;
    double H, S, I, a, b, c;

    for (i = 0; i < nl * nc; i++) {

        a = (double) R[i] / mn;
        b = (double) G[i] / mn;
        c = (double) B[i] / mn;

        rbg2hsi(a, b, c, &H, &S, &I);

        // Zerar H e S para deixar em tons de cinza
        H = 0;
        S = 0;

        //transforma HSI >> RGB
        hsi2rgb(H, S, I, &a, &b, &c);
        R[i] = a * mn;
        G[i] = b * mn;
        B[i] = c * mn;

    }
}

void msg(char *s) {
    printf("\nWatershed");
    printf("\n---------");
    printf("\nUSO.:  %s  <IMG.PGM>", s);
    printf("\nONDE:\n");
    printf("    <IMG.PGM>  Arquivo da imagem em formato PGM ASCII\n\n");
    exit(10);
}

/*------------------------------------------------------
 *          P R O G R A M A    P R I N C I P A L
 *------------------------------------------------------*/
int main(int argc, char *argv[]) {
    int nc, nl, mn;
    char nomepgm[100];
    char nomeppm[100];
    char comando[100];

    image R, G, B, Rc, Gc, Bc;
    image Out, Grd;
    if (argc < 2)
        msg(argv[0]);

    read_ppm(&R, &G, &B, argv[1], &nl, &nc, &mn);
    img_info(argv[1], nl, nc, mn);
    Rc = img_alloc(nl, nc);
    Gc = img_alloc(nl, nc);
    Bc = img_alloc(nl, nc);
    
    for(int i = 0; i< nl * nc; i++){
        Rc[i] = R[i];
        Gc[i] = G[i];
        Bc[i] = B[i];
    }
    
    Out = img_alloc(nl, nc);
    Grd = img_alloc(nl, nc);

    tonsCinza(Rc, Gc, Bc, nl, nc, mn);
    
    gradient(Rc, Grd, nl, nc, mn, 1);
    watershed(Grd, Out, nl, nc, mn, (nl / 2), nc / 2);

    tonsCinzaParc(R, G, B, Out, nl, nc, mn);

    sprintf(nomepgm, "%s-%s", argv[1], "watershed.pgm");
    write_pgm(Out, nomepgm, nl, nc, mn);

    sprintf(nomeppm, "%s-%s", argv[1], "cor-pb.ppm");
    write_ppm(R, G, B, nomeppm, nl, nc, mn);

    strcpy(comando, "./i_view32 "); // comando pra executar no linux
    strcat(comando, nomeppm);
    strcat(comando, " /hide=7");
    system(comando);

    img_free(R);
    img_free(G);
    img_free(B);
    img_free(Rc);
    img_free(Gc);
    img_free(Bc);
    img_free(Out);
    img_free(Grd);
    return 0;
}