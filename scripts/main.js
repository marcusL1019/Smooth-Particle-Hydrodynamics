//-------------------------Table Of Contents-----------------------------
//---VectorClass
//---ArrayFunctions
//---ParticleClass
//---CellClass
//---GridClass
//---SPHClass
//------constructor
//------wPolyFunction
//------wPolyGradientFunction
//------wPolyLaplacianFunction
//------wPressGradientFunction
//------wViscLaplacianFunction
//------updateNeighbor
//------updateRho
//------updateForces
//------integrate
//------boundary
//------draw
//------update

//-------------------------VectorClass----------------------------
//vector class to simplify vector math
class Vector {
	constructor(x, y) {
		this.x = x;
		this.y = y;
	}

	set(a, b) {
		this.y = a;
		this.y = b;
	}

	clone() {
		return new Vector(this.x, this.y);
	}

	scalarAdd(a, b) {
		this.set(this.x + a, this.y + b);
	}

	scalarMult(a, b) {
		this.x *= a;
		this.y *= b;
	}

	vectorAdd(v) {
		this.x += v.x;
		this.y += v.y;
	}

	vectorSub(v) {
		this.x -= v.x;
		this.y -= v.y;
	}

	vectorMult(v) {
		this.x *= v.x;
		this.y *= v.y;
	}

	distance(v) {
		return Math.sqrt(
			(this.x - v.x) * (this.x - v.x) + (this.y - v.y) * (this.y - v.y)
		);
	}

	rSqr(v) {
		return (this.x - v.x) * (this.x - v.x) + (this.y - v.y) * (this.y - v.y);
	}

	unit() {
		var r = Math.sqrt(this.x * this.x + this.y * this.y);
		var x = this.x / r;
		var y = this.y / r;
		return new Vector(x, y);
	}
}

//-------------------------ArraysFunctions----------------------------
//functions for creating 1D/2D arrays
var array1D = function (n) {
	var out = new Array();
	out.length = n;

	return out;
};

var array2D = function (n, m) {
	var N = new Array();
	N.length = n;

	for (var i = 0; i < n; i++) {
		var M = new Array();
		M.length = m;
		M.fill(0);
		N[i] = M;
	}

	return N;
};

//-------------------------ParticleClass----------------------------

//particle object stores id, position, velocity, force, pressure, and density
//of particle.
class Particle {
	constructor(id, x, y) {
		this.id = id;
		this.pos = new Vector(x, y);
		this.v = new Vector(0, 0);
		this.f = new Vector(0, 0);
		this.rho = 0;
		this.p = 0;
	}
}

//-------------------------CellClass----------------------------
class Cell {
	constructor(capacity) {
		if (capacity == undefined) capacity = 0;
		this.values = array1D(capacity);
		this.count = capacity;
	}

	clear() {
		this.values.length = 0;
		this.count = 0;
	}

	add(id) {
		this.values.push(id);
		this.count++;
	}
}

//-------------------------GridClass----------------------------

class Grid {
	constructor(width, height, numDivisions) {
		this.width = width;
		this.height = height;
		this.numDivisions = numDivisions;

		this.gridX = this.width / this.numDivisions;
		this.gridY = this.height / this.numDivisions;

		this.cells = array2D(this.numDivisions, this.numDivisions);

		for (var i = 0; i < this.numDivisions; i++) {
			for (var j = 0; j < this.numDivisions; j++) {
				this.cells[i][j] = new Cell(0);
			}
		}
	}

	clear() {
		for (var i = 0; i < this.numDivisions; i++) {
			for (var j = 0; j < this.numDivisions; j++) {
				this.cells[i][j].clear();
			}
		}
	}

	add(position, id) {
		var x = Math.floor(position.x / this.gridX);
		var y = Math.floor(position.y / this.gridY);

		if (x < 0) x = 0;
		if (y < 0) y = 0;
		if (x >= this.numDivisions) x = this.numDivisions - 1;
		if (y >= this.numDivisions) y = this.numDivisions - 1;

		this.cells[y][x].add(id);
	}
}

//-------------------------SPHClass----------------------------

var dt = 0.5;
var canvas = document.getElementById('particleCanvas');
var numDivisions = 100;
var material = {
	m: 0, //particle mass
	k: 15, //gas constant
	mu: 150, //dynamic viscosity
	rho0: 1000, //rest density
	sigma: 10, //surface tension parameter
	damping: 0.1, //boundary condition velocity dampening
};

class SPH {
	constructor(canvas, numDivisions, material, dt) {
		this.canvas = canvas;
		this.context = canvas.getContext('2d');
		this.width = canvas.width;
		this.height = canvas.height;
		this.gridSizeX = this.width / numDivisions;
		this.gridSizeY = this.height / numDivisions;
		this.material = material;
		this.dt = dt;
		this.numDivisions = numDivisions;
		this.damping = this.material.damping;
		this.maxRho = 0;

		//compute simulation parameters
		this.area = this.width * this.height;
		this.numParticles = Math.floor(
			Math.pow((numDivisions + numDivisions) / 4, 2)
		);
		this.material.m = (this.material.rho0 * this.area) / this.numParticles;
		this.h = Math.min(this.gridSizeX, this.gridSizeY);

		//generate particles
		this.particles = array1D(this.numParticles);
		this.grid = new Grid(this.width, this.height, numDivisions);

		var x = this.h;
		//var x = this.width / 4;
		var y = 25 * this.h;
		var del = this.h;

		for (var i = 0; i < this.numParticles; i++) {
			if (y > this.height - 2 * this.h) {
				y = 25 * this.h;
				x += del;
			}
			this.particles[i] = new Particle(i, x + Math.random(), y + Math.random());

			y += del;
		}

		//precompute kernel scalar components and other useful
		//quantities
		this.h2 = this.h * this.h;
		this.wPoly = 315 / (64 * Math.PI * Math.pow(this.h, 9));
		this.wPolyGradient = -945 / (32 * Math.PI * Math.pow(this.h, 9));
		this.wPolyLaplacian = this.wPoly;
		this.wPressGradient = -45 / (Math.PI * Math.pow(this.h, 6));
		this.wViscLaplacian = -1 * this.wPressGradient;

		//call update method to begin animation
		this.update();
	}

	wPolyFun(r2) {
		var dif = Math.abs(this.h2 - r2);
		return this.wPoly * Math.pow(dif, 3);
	}

	wPolyGradientFun(r2, v) {
		var dif = (this.h2 - r2) * (this.h2 - r2);
		dif *= this.wPolyGradient;
		v.scalarMult(dif, dif);
		return new Vector(v.x, v.y);
	}

	wPolyLaplacianFun(r2) {
		var dif = (this.h2 - r2) * (7 * r2 - 3 * this.h2);
		return this.wPolyLaplacian * dif;
	}

	wPressGradientFun(r, vUnit) {
		var dif = (this.h - r) * (this.h - r);
		dif *= this.wPressGradient;
		vUnit.scalarMult(dif, dif);
		return new Vector(vUnit.x, vUnit.y);
	}

	wViscLaplacianFun(r) {
		var dif = this.h - r;
		return dif * this.wViscLaplacian;
	}

	updateNeighbor() {
		//console.log('neighbor');

		this.grid.clear(); //clear grid before adding particles to current cells

		for (var i = 0; i < this.numParticles; i++) {
			var pi = this.particles[i];
			//if (i == 28) console.log(pi.pos.x, pi.pos.y);
			this.grid.add(pi.pos, pi.id); //add particles to the correct cells

			//set forces, p, and density to zero
			pi.f.set(0, 0);
			pi.p = 0;
			pi.rho = 0;
		}
	}

	updateRho() {
		//console.log('rho');
		this.maxRho = 0;
		for (var i = 0; i < this.numParticles; i++) {
			pi = this.particles[i];

			//compute particles' density contribution from itself
			pi.rho += this.material.m * this.wPolyFun(0);

			//compute bounds of particle (pos +/- h) to determine which cells to look at for neighbors
			var minX = Math.floor((pi.pos.x - this.h) / this.gridSizeX);
			var minY = Math.floor((pi.pos.y - this.h) / this.gridSizeY);
			var maxX = Math.floor((pi.pos.x + this.h) / this.gridSizeX);
			var maxY = Math.floor((pi.pos.y + this.h) / this.gridSizeY);

			//check if particle closer to edge than h
			if (minX < 0) minX = 0;
			if (minY < 0) minY = 0;
			if (maxX >= this.numDivisions) maxX = this.numDivisions - 1;
			if (maxY >= this.numDivisions) maxY = this.numDivisions - 1;

			//loop through cells that might contain neighbors
			for (var y = minY; y <= maxY; y++) {
				for (var x = minX; x <= maxX; x++) {
					if (this.grid.cells[y] == undefined)
						console.log(i, this.numDivisions, maxY);
					var cell = this.grid.cells[y][x];

					for (var j = 0; j < cell.count; j++) {
						var id_j = cell.values[j];

						//if pj != pi, compute distance and density of pi
						if (id_j != i) {
							var pj = this.particles[id_j];

							var r2 = pi.pos.rSqr(pj.pos);

							if (r2 <= this.h2) {
								var w = this.wPolyFun(r2);
								pi.rho += this.material.m * w;
								if (pi.rho > this.maxRho) this.maxRho = pi.rho;
							}
						}
					}
				}
			}
		}

		//update pressures
		for (var i = 0; i < this.numParticles; i++) {
			var pi = this.particles[i];

			pi.p = this.material.k * (pi.rho - this.material.rho0);
		}
	}

	updateForces() {
		//console.log('forces');

		for (var i = 0; i < this.numParticles; i++) {
			var pi = this.particles[i];
			pi.f.set(0, 0.1 * 9.8 * pi.rho);

			var color = new Vector(0, 0);
			var colorLaplacian = 0;

			//compute bounds of particle (pos +/- h) to determine which cells to look at for neighbors
			var minX = Math.floor((pi.pos.x - this.h) / this.gridSizeX);
			var minY = Math.floor((pi.pos.y - this.h) / this.gridSizeY);
			var maxX = Math.floor((pi.pos.x + this.h) / this.gridSizeX);
			var maxY = Math.floor((pi.pos.y + this.h) / this.gridSizeY);

			//check if particle closer to edge than h
			if (minX < 0) minX = 0;
			if (minY < 0) minY = 0;
			if (maxX >= this.numDivisions) maxX = this.numDivisions - 1;
			if (maxY >= this.numDivisions) maxY = this.numDivisions - 1;

			//loop through cells that might contain neighbors
			for (var y = minY; y <= maxY; y++) {
				for (var x = minX; x <= maxX; x++) {
					var cell = this.grid.cells[y][x];

					for (var j = 0; j < cell.count; j++) {
						var id_j = cell.values[j];

						//if pj != pi, compute distance and density of pi
						if (id_j != i) {
							var pj = this.particles[id_j];

							var r = pi.pos.distance(pj.pos);
							var r2 = r * r;

							if (r < this.h) {
								//compute pressure forces
								var vUnit = pi.pos.clone();
								vUnit.vectorSub(pj.pos);
								vUnit = vUnit.unit();
								var wPress = this.wPressGradientFun(r, vUnit);
								var a =
									(this.material.m / pi.rho) *
									(this.material.m / pj.rho) *
									((pi.p + pj.p) / 2);
								wPress.scalarMult(a, a);
								pi.f.vectorSub(wPress);

								//compute viscosity forces
								var wVisc = this.wViscLaplacianFun(r);
								var vVisc = pj.v.clone();
								vVisc.vectorSub(pi.v);
								var b =
									(this.material.m / pi.rho) *
									(this.material.m / pj.rho) *
									this.material.mu *
									wVisc;
								vVisc.scalarMult(b, b);
								pi.f.vectorAdd(vVisc);

								//compute color field gradient
								var vColor = pi.pos.clone();
								vColor.vectorSub(pj.pos);
								vColor = this.wPolyGradientFun(r2, vColor);
								var v = this.material.m / pj.rho;
								vColor.scalarMult(v, v);
								color.vectorAdd(vColor);
								color = color.unit();

								//compute color field laplacian
								var wC2 = this.wPolyLaplacianFun(r2) * v;
								colorLaplacian += wC2;

								//compute surface tension
								color.scalarMult(
									colorLaplacian * this.material.sigma,
									colorLaplacian * this.material.sigma
								);
							}
						}
					}
				}
			}
			pi.f.vectorSub(color);
		}
	}

	integrate() {
		//console.log('integrate');

		for (var i = 0; i < this.numParticles; i++) {
			var pi = this.particles[i];
			//if (i == 200) console.log(pi.v.x, pi.v.y);
			pi.f.scalarMult(this.dt / this.material.m, this.dt / this.material.m);

			pi.v.vectorAdd(pi.f);
			pi.pos.vectorAdd(pi.v);
		}
	}

	boundary() {
		for (var i = 0; i < this.numParticles; i++) {
			var pi = this.particles[i];

			if (pi.pos.y > this.canvas.height - this.h) {
				pi.v.y = -this.damping * Math.abs(pi.v.y);
				pi.pos.y = this.canvas.height - this.h;
			}
			if (pi.pos.y < this.h) {
				pi.v.y = this.damping * Math.abs(pi.v.y);
				pi.pos.y = this.h;
			}
			if (pi.pos.x > this.canvas.width - this.h) {
				pi.v.x = -this.damping * Math.abs(pi.v.x);
				pi.pos.x = this.canvas.width - this.h;
			}
			if (pi.pos.x < this.h) {
				pi.v.x = this.damping * Math.abs(pi.v.x);
				pi.pos.x = this.h;
			}
		}
	}

	draw() {
		//console.log('draw');

		for (var i = 0; i < this.numParticles; i++) {
			var pi = this.particles[i];
			var transparency = pi.rho / this.maxRho;

			this.context.beginPath();
			this.context.arc(
				this.particles[i].pos.x,
				this.particles[i].pos.y,
				this.h / 3,
				0,
				2 * Math.PI
			);
			this.context.fillStyle = 'rgba(175, 238, 238,' + transparency + ')';
			this.context.fill();
		}
	}
	update() {
		this.context.clearRect(0, 0, this.width, this.height);

		this.draw();
		this.updateNeighbor();
		this.updateRho();
		this.updateForces();
		this.integrate();
		this.boundary();

		requestAnimationFrame(this.update.bind(this));
	}
}

var sph = new SPH(canvas, numDivisions, material, dt);
