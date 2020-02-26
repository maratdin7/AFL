#include <iostream>
#include <unistd.h>
#include <fcntl.h>

using namespace std;

int main() {
    char[2] a;
        const char *fn_a = "example/simple_test/test_a";
        int fd_a = open(fn, O_CREAT | O_RDWR, 0600);
	a = read(fd_a, a, 2); 
    int b = 1;

    if (a == "10") {
        const char *fn = "example/simple_test/test";
        int fd = open(fn, O_CREAT | O_RDWR, 0600);
        if (fd >= 0) {
            char ve[2];

            if (read(fd, ve, 2) == 2) {
                while (__AFL_LOOP(1000)) {
                    if (ve[0] == 'v') {
                        a = 15;
                        if (ve[1] == 'e') b--;
                    }
                }
            }
        }
    }
    cout << a / b;
    return 0;
}
