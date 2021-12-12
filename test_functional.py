from MyBWT import MyBWT
from ag2 import BWTstore
from bwtclass import BWTtrans
import generator
import ag_tra


# def test_generate_LC():
#     for i in range(2, 5):
#         for j in range(5):
#             ref = generator.generate_DNA(10 ** i)
#             my_lc = "".join(MyBWT(ref).generate_LC(ref))
#             lc = ag_tra.fl(ref)
#             print(my_lc == lc)
#             assert my_lc == lc


# if __name__ == "__main__":
#     test_generate_LC()

# print(BWTtrans("123", "2").generateTally())

r = generator.generate_DNA(100)
sr = generator.extract(r, 10, 1)
# print(MyBWT("ASDFG").query("SD"))
print(MyBWT(r).query(sr[0][0]), sr[1][0])
# print(sr[1])
