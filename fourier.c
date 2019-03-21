#include <math.h>

#include "fourier.h"

void nft(double complex s[MAX_SIZE], double complex t[MAX_SIZE], int n, int sign) {
  for (int k = 0; k < n; k++) {
    t[k] = 0;

    for (int j = 0; j < n; j++) {
      t[k] += s[j] * cexp(sign * 2 * PI * k * j * I / n);
    }
  }
}

void nft_forward(double complex s[MAX_SIZE], double complex t[MAX_SIZE], int n) {
  nft(s, t, n, -1);
}

void nft_inverse(double complex t[MAX_SIZE], double complex s[MAX_SIZE], int n) {
  nft(t, s, n, 1);

  for (int k = 0; k < n; k++) {
    s[k] /= n;
  }
}

void fft(double complex s[MAX_SIZE], double complex t[MAX_SIZE], int n, int sign) {

  double complex pares[n / 2];
  double complex impares[n / 2];
  double complex paresT[n / 2];
  double complex imparesT[n / 2];

  //Caso base
  if (n == 1) {
    t[0] = s[0];
    return;
  }

  //Dividir a entrada original em partes menores (separar em impares e pares para a recursao)
  int a = 0, b = 0;
  for (int i = 0; i < n; i++) {
    //Se retornar 1 eh impar (resto da divisao) e ja serve como True
    if (i % 2) {
      impares[a++] = s[i];
    } else {
      pares[b++] = s[i];
    }
  }

  //Resolver recursivamente o problema para cada parte
  fft(pares, paresT, n / 2, sign);
  fft(impares, imparesT, n / 2, sign);

  //Combinar as solucoes para resolver o problema para a entrada original
  for (int i = 0; i < n / 2; i++) {
    t[i] = paresT[i] + imparesT[i] * cexp(sign * 2 * PI * i * I / n);
    t[i + n / 2] = paresT[i] - imparesT[i] * cexp(sign * 2 * PI * i * I / n);
  }
}

void fft_forward(double complex s[MAX_SIZE], double complex t[MAX_SIZE], int n) {
  fft(s, t, n, -1);
}

void fft_inverse(double complex t[MAX_SIZE], double complex s[MAX_SIZE], int n) {
  fft(t, s, n, 1);

  for (int k = 0; k < n; k++) {
    s[k] /= n;
  }
}

void fft_forward_2d(double complex matrix[MAX_SIZE][MAX_SIZE], int width, int height) {
  // Cria os vetores das linhas e colunas de T e S
  double complex tline[width];
  double complex sline[width];
  double complex tcolumn[height];
  double complex scolumn[height];

  for (int i = 0; i < height; i++) { // para cada linha da matriz
    for (int j = 0; j < width; j++) { // para cada elemento da linha
      sline[j] = matrix[i][j]; // cada linha tem um sline com os valores da matriz original
    }
    fft_forward(sline, tline, width); // aplique a transformada unidimensional normal sobre cada linha
    for (int j = 0; j < width; j++) {
      matrix[i][j] = tline[j]; // guarda na matriz original
    }
  }
  for (int j = 0; j < width; j++) { // para cada coluna da matriz
    for (int i = 0; i < height; i++) {  // para cada elemento da coluna
      scolumn[i] = matrix[i][j]; // cada coluna tem um scolumn com os valores da matriz original
    }
    fft_forward(scolumn, tcolumn, height); // aplique a transformada unidimensional normal sobre a coluna
    for (int i = 0; i < height; i++) {
      matrix[i][j] = tcolumn[i]; // guarda na matriz original
    }
  }
  return;
}

void fft_inverse_2d(double complex matrix[MAX_SIZE][MAX_SIZE], int width, int height) {
  // Cria os vetores das linhas e colunas de T e S
  double complex tline[width];
  double complex sline[width];
  double complex tcolumn[height];
  double complex scolumn[height];

  for (int i = 0; i < height; i++) { // para cada linha da matriz
    for (int j = 0; j < width; j++) { // para cada elemento da linha
      sline[j] = matrix[i][j]; // cada linha tem um sline com os valores da matriz original
    }
    fft_inverse(sline, tline, width); // aplique a transformada unidimensional inversa normal sobre cada linha
    for (int j = 0; j < width; j++) {
      matrix[i][j] = tline[j]; // guarda na matriz original
    }
  }
  for (int j = 0; j < width; j++) { // para cada coluna da matriz
    for (int i = 0; i < height; i++) { // para cada elemento da coluna
      scolumn[i] = matrix[i][j]; // cada coluna tem um scolumn com os valores da matriz original
    }
    fft_inverse(scolumn, tcolumn, height); // aplique a transformada unidimensional inversa normal sobre a coluna
    for (int i = 0; i < height; i++) {
      matrix[i][j] = tcolumn[i]; // guarda na matriz original
    }
  }
  return;
}

void filter(double complex input[MAX_SIZE][MAX_SIZE], double complex output[MAX_SIZE][MAX_SIZE], int width, int height, int flip) {
  int center_x = width / 2;
  int center_y = height / 2;

  double variance = -2 * SIGMA * SIGMA;

  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      int dx = center_x - (x + center_x) % width;
      int dy = center_y - (y + center_y) % height;

      double d = dx * dx + dy * dy;

      double g = exp(d / variance);

      if (flip) {
        g = 1 - g;
      }

      output[y][x] = g * input[y][x];
    }
  }
}

void filter_lp(double complex input[MAX_SIZE][MAX_SIZE], double complex output[MAX_SIZE][MAX_SIZE], int width, int height) {
  filter(input, output, width, height, 0);
}

void filter_hp(double complex input[MAX_SIZE][MAX_SIZE], double complex output[MAX_SIZE][MAX_SIZE], int width, int height) {
  filter(input, output, width, height, 1);
}
