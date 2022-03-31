#include <tr1/memory>

struct Foo
{
    int a;
};

std::tr1::shared_ptr<Foo> p;

extern "C" int plus_one(int a) {
    p.reset(new Foo);
    return a + 1;
}
