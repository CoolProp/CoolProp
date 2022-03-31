extern "C" void myfunc(int a, double b, double c[]);

int main() {
    double c[20];

    for (int i = 0; i < 20; i++)
        c[i] = double(i);

    myfunc(1, 42.0, c);
}