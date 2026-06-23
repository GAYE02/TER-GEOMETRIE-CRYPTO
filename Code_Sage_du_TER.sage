# Decompoisition de Mumford
def mumford_decomposition(P1,P2):
				x1, y1 = P1[0], P1[1]
				x2, y2 = P2[0], P2[1]
				u = R([x1*x2,-(x1+x2), 1])
				u = u.monic()
				v = R([(y1-y2)/(x1-x2)*x1 + y1, (y1-y2)/(x1-x2)])
				if (v^2 - f) // u == 0:
					print("v^2 - f n'est pas multiple de u")
				return (u, v)

# Addition de Cantor
def addition_Cantor(D1,D2,f):
			u1,v1 = D1[0],D1[1]
			u2,v2 = D2[0],D2[1]
			d1,e1,e2 =  u1.xgcd(u2)
			d,c1,c2 =  d1.xgcd(v1+v2)
			s1 = c1*e1
			s2 = c1*e2
			s3 = c2
			u = (u1*u2) // pow(d,2) 
			v = (s1*u1*v2 + s2*u2*v1 + s3(v1*v2 + f)) % u
			return (u,v)
# Reduction de Cantor
def reduction_Cantor(D, f):
				u, v = D[0], D[1]
				g = (f.degree() - 1) // 2
				while u.degree() > g:
						U = (f - v^2) // u
						v = -(v) % U
						u = U
						v = V
				return (U, V)

# Decomposition de polynomial (Richelot)
def decomposition_polynomial(G1,G2,G3):
				G1_list = G1.list()
				G2_list = G2.list()
				G3_list = G3.list()
				Matrice = Matrix([G1_list,G2_list,G3_list])
				Delta = det(Matrice)
				if Delta != 0:
						H1 = G2.derivative(x)*G3 - G2*G3.derivative(x)
						H2 = G3.derivative(x)*G1 - G3*G1.derivative(x)
						H3 = G1.derivative(x)*G2 - G1*G2.derivative(x)
						return Delta^-1 * R(H1) * R(H2) * R(H3)
				else: return "Delta est nul"


