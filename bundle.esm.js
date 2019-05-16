var Vector2 = function(x, y) {
		if(typeof y === "undefined")
			throw Error("initialize Vector2 with less than 2 parameter");

		this.x = x;
		this.y = y;
	};

	Vector2.Zero = new Vector2(0, 0);
	Vector2.One = new Vector2(1, 1);
	Vector2.XAxis = new Vector2(1, 0);
	Vector2.YAxis = new Vector2(0, 1);
	Vector2.fromAngle = function(theta) {
		return new Vector2(Math.cos(theta), Math.sin(theta));
	};
	
	Vector2.prototype.copy = function() {
		return new Vector2(this.x, this.y);
	};

	Vector2.prototype.set = function(vector) {
		this.x = vector.x;
		this.y = vector.y;
		return this;
	};

	Vector2.prototype.setXY = function(x, y) {
		this.x = x;
		this.y = y;
		return this;
	};

	Vector2.prototype.clear = function() {
		this.x = 0;
		this.y = 0;
	};

	Vector2.prototype.floor = function() {
		return new Vector2(
			Math.floor(this.x),
			Math.floor(this.y)
		);
	};

	Vector2.prototype.floorSelf = function() {
		this.x = Math.floor(this.x);
		this.y = Math.floor(this.y);
	};

	Vector2.prototype.frac = function() {
		return new Vector2(
			this.x - Math.floor(this.x),
			this.y - Math.floor(this.y)
		);
	};

	Vector2.prototype.fracSelf = function() {
		this.x -= Math.floor(this.x);
		this.y -= Math.floor(this.y);
	};

	Vector2.prototype.equals = function(vector) {
		return this.x == vector.x && this.y == vector.y;
	};

	Vector2.prototype.notEquals = function(vector) {
		return this.x != vector.x || this.y != vector.y;
	};

	Vector2.prototype.add = function(vector) {
		return new Vector2(
			this.x + vector.x,
			this.y + vector.y
		);
	};
	Vector2.prototype.addSelf = function(vector) {
		this.x += vector.x;
		this.y += vector.y;
		return this;
	};
		
	Vector2.prototype.weightedAddSelf = function(vector, scalar) {
		this.x += vector.x * scalar;
		this.y += vector.y * scalar;
		return this;
	};
		
	Vector2.prototype.sub = function(vector) {
		return new Vector2(
			this.x - vector.x,
			this.y - vector.y
		);
	};
	Vector2.prototype.subSelf = function(vector) {
		this.x -= vector.x;
		this.y -= vector.y;
		return this;
	};
		
	// scaling!
	Vector2.prototype.mulFloat = function(right) {
		return new Vector2(
			this.x * right,
			this.y * right
		);
	};
	Vector2.prototype.mulFloatSelf = function(right) {
		this.x *= right;
		this.y *= right;
		return this;
	};
		
	Vector2.prototype.divFloat = function(right) {
		var inv = 1.0 / right;
		return new Vector2(
			this.x * inv,
			this.y * inv
		);
	};
	Vector2.prototype.divFloatSelf = function(right) {
		this.x /= right;
		this.y /= right;
		return this;
	};

	// per-element multiplication
	Vector2.prototype.mulVector = function(right) {
		return new Vector2(
			this.x * right.x,
			this.y * right.y
		);
	};
	Vector2.prototype.mulVectorSelf = function(right) {
		this.x *= right.x;
		this.y *= right.y;
		return this;
	};

	Vector2.prototype.divVector = function(right) {
		return new Vector2(
			this.x / right.x,
			this.y / right.y
		);
	};
	Vector2.prototype.divVectorSelf = function(right) {
		this.x /= right.x;
		this.y /= right.y;
		return this;
	};

	Vector2.prototype.positive = function() { return this; };
	Vector2.prototype.negative = function() {
		return new Vector2(-this.x, -this.y);
	};

	// helpers

	Vector2.prototype.length = function() {
		return Math.sqrt(this.x*this.x + this.y*this.y);
	};
	Vector2.prototype.lengthSquared = function() {
		return this.x*this.x + this.y*this.y;
	};
	Vector2.prototype.distance = function(right) {
		var x = this.x - right.x;
		var y = this.y - right.y;
		return Math.sqrt(x*x + y*y);
	};
	Vector2.prototype.distanceSquared = function(right) {
		var x = this.x - right.x;
		var y = this.y - right.y;
		return x*x + y*y;
	};
	Vector2.prototype.normalize = function() {
		var length = Math.sqrt(this.x*this.x + this.y*this.y);
		if(length > 1e-08) {
			var invL = 1.0 / length;
			this.x *= invL;
			this.y *= invL;
		}
		return length;
	};

	Vector2.prototype.normalizedCopy = function() {
		var ret = this.copy();
		ret.normalize();
		return ret;
	};

	Vector2.prototype.dotProduct = function(vector) {
		return this.x*vector.x + this.y*vector.y;
	};

	Vector2.prototype.getPerpendicular = function() {
		return this.getRightPerpendicular();
	};

	Vector2.prototype.getLeftPerpendicular = function() {
		var x = this.y;
		var y = -1 * this.x;
		return new Vector2(x, y);
	};

	Vector2.prototype.getRightPerpendicular = function() {
		var x = -1 * this.y;
		var y = this.x;
		return new Vector2(x, y);
	};

	Vector2.prototype.makePerpendicular = function() {
		var tempX = this.x;
		this.x = -this.y;
		this.y = tempX;
	};

	Vector2.prototype.crossProduct = function(vector) {
		return this.x * vector.y + this.y * vector.x;
	};

	Vector2.prototype.lerp = function(to, i) {
		return this.add(to.sub(this).mulFloat(i));
	};

	Vector2.prototype.lerpSelf = function(to, i) {
		return this.weightedAddSelf(to.sub(this), i);
	};

	Vector2.prototype.slerp = function(to, i) {
		return this.add(to.sub(this).mulFloat(
			0.5 + (Math.sin((3.141592654 * i) - 1.570796) * 0.5)));
	};

	Vector2.prototype.rotate = function(theta) {
		var co = Math.cos(theta);
		var si = Math.sin(theta);
		return new Vector2(
			co * this.x - si * this.y,
			si * this.x + co * this.y
		);
	};

	Vector2.prototype.rotateSelf = function(theta) {
		var co = Math.cos(theta);
		var si = Math.sin(theta);
		var xx = co * this.x - si * this.y;
		this.y = si * this.x + co * this.y;
		this.x = xx;
	};
	
	// get (signed and directional) angle between this and the given vector in degrees 
	Vector2.prototype.getDirectedAngle = function(point) {
		return Math.atan2(this.crossProduct(point), this.dotProduct(point)) * 180 / Math.PI;
	};

	Vector2.prototype.reflectOnNormal = function(normal) {
		//v' = 2 * (v . n) * n - v
		var newVector =
			this.sub(
				normal
				.mulFloat(this.dotProduct(normal))
				.mulFloat(2)
			);
		return newVector;
		
	};

	Vector2.prototype.toCartesian = function() {
		return new Vector2(
			this.x * Math.cos(this.y),
			this.x * Math.sin(this.y)
		);
	};

	Vector2.prototype.toPolar = function() {
		return new Vector2(
			Math.sqrt(this.x * this.x + this.y * this.y),
			Math.atan2(this.y, this.x)
		);
	};

	Vector2.prototype.signum = function() {
		return new Vector2(
			this.x.sign(),
			this.y.sign()
		);
	};

	Vector2.prototype.absolute = function() {
		return new Vector2(
			Math.abs(this.x),
			Math.abs(this.y)
		);
	};


	Vector2.prototype.toJson = function() {
		var resultJson = {
			"x": this.x,
			"y": this.y
		};

		return resultJson;
	};
	
	Vector2.fromJson = function(vectorJson) {
		return new Vector2(vectorJson.x, vectorJson.y);
	};

var Invalid = true;
	var Valid = false;

	var AABB = function(minPt, maxPt) {
		this.Min = minPt || Vector2.Zero.copy();
		this.Max = maxPt || Vector2.Zero.copy();
		this.Validity = typeof minPt === "undefined" ? Invalid : Valid;
	};	

	AABB.prototype.clear = function() {
		this.Min.set(Vector2.Zero);
		this.Max.set(Vector2.Zero);
		this.Validity = Invalid;
	};	

	AABB.prototype.expandToInclude = function(pt) {
		if (this.Validity == Valid) {
			if (pt.x < this.Min.x) {
				this.Min.x = pt.x;
			} else if (pt.x > this.Max.x) {
				this.Max.x = pt.x;
			}
			
			if (pt.y < this.Min.y) {
				this.Min.y = pt.y;
			} else if (pt.y > this.Max.y) {
				this.Max.y = pt.y;
			}
		} else {
			this.Min.set(pt);
			this.Max.set(pt);
			this.Validity = Valid;
		}	};	

	AABB.prototype.contains = function(pt) {
		if (this.Validity == Invalid) { return false; }
		return ((pt.x >= this.Min.x) && (pt.x <= this.Max.x) && (pt.y >= this.Min.y) && (pt.y <= this.Max.y));
	};
	
	AABB.prototype.containsAABB = function(other) {
		if (this.Validity == Invalid) { return false; }		if (other.Validity == Invalid) { return false; }		
		return (other.Min.x >= this.Min.x) &&
			(other.Max.x <= this.Max.x) &&
			(other.Min.y >= this.Min.y) &&
			(other.Max.y <= this.Max.y);
	};
	
	AABB.prototype.intersects = function(box) {
			var overlapX = ((this.Min.x <= box.Max.x) && (this.Max.x >= box.Min.x));
			var overlapY = ((this.Min.y <= box.Max.y) && (this.Max.y >= box.Min.y));
			
			return (overlapX && overlapY);
	};	

	AABB.prototype.getSize = function() { return this.Max.sub(this.Min); };

	AABB.prototype.getMiddle = function() {
		return this.Min.add(this.Max.sub(this.Min).mulFloat(0.5));
	};
	AABB.prototype.getTopLeft = function() {
		return this.Min;
	};
	AABB.prototype.getTopRight = function() {
		return new Vector2(this.Max.x, this.Min.y);
	};
	AABB.prototype.getBottomLeft = function() {
		return new Vector2(this.Min.x, this.Max.y);
	};
	AABB.prototype.getBottomRight = function() {
		return this.Max;
	};
	AABB.prototype.subdivide = function() {
		var min = this.Min,
			middle = this.getMiddle(),
			max = this.Max;
		
		var i = new AABB(min, middle);
		return [
		    new AABB(min, middle),
		    new AABB(new Vector2(middle.x, min.y), new Vector2(max.x, middle.y)),
		    new AABB(new Vector2(min.x, middle.y), new Vector2(middle.x, max.y)),
		    new AABB(middle, max)
		];
	};

var Material = function(materialIndex) {
  this.setColor(Material.getColor(materialIndex));

  this.materialIndex = materialIndex;
  this.particleMass = 1;
  this.restDensity = 1;
  this.stiffness = 1;
  this.bulkViscosity = 1;
  this.surfaceTension = 0.2;
  this.elasticity = 0.05;
  this.maxDeformation = 0;
  this.meltRate = 0;
  this.shearViscosity = 0;
  this.viscosity = 0.02;
  this.damping = 0.1;
  this.wallFriction = 0;
  this.wallAttraction = 1;
  this.smoothing = 1;
  this.isElastic = false;
  this.springK = 0.3;

  this.yieldPoint = 0;
  this.yieldRate = 1;
};

// debug colors
Material.getColor = function(index) {
  var materialColors = [
    "#1f78b4",
    "#e31a1c",
    "#fdbf6f",
    "#b2df8a",
    "#fb9a99",
    "#ff7f00",
    "#33a02c",
    "#cab2d6",
    "#6a3d9a",
    "#ffff99",
    "#b15928"
  ];

  return materialColors[index % materialColors.length];
};

// Property setters
Material.prototype.setColor = function(color) {
  this.color = color;

  return this;
};

Material.prototype.setParticleMass = function(particleMass) {
  this.particleMass = particleMass;

  return this;
};

Material.prototype.setRestDensity = function(restDensity) {
  this.restDensity = restDensity;

  return this;
};

Material.prototype.setStiffness = function(stiffness) {
  this.stiffness = stiffness;

  return this;
};

Material.prototype.setBulkViscosity = function(bulkViscosity) {
  this.bulkViscosity = bulkViscosity;

  return this;
};

Material.prototype.setSurfaceTension = function(surfaceTension) {
  this.surfaceTension = surfaceTension;

  return this;
};

Material.prototype.setElasticity = function(elasticity) {
  this.elasticity = elasticity;

  return this;
};

Material.prototype.setMaxDeformation = function(maxDeformation) {
  this.maxDeformation = maxDeformation;

  return this;
};

Material.prototype.setMeltRate = function(meltRate) {
  this.meltRate = meltRate;

  return this;
};

Material.prototype.setShearViscosity = function(shearViscosity) {
  this.shearViscosity = shearViscosity;

  return this;
};

Material.prototype.setViscosity = function(viscosity) {
  this.viscosity = viscosity;

  return this;
};

Material.prototype.setDamping = function(damping) {
  this.damping = damping;

  return this;
};

Material.prototype.setWallFrictiong = function(wallFriction) {
  this.wallFriction = wallFriction;

  return this;
};

Material.prototype.setWallAttraction = function(wallAttraction) {
  this.wallAttraction = wallAttraction;

  return this;
};

Material.prototype.setSmoothing = function(smoothing) {
  this.smoothing = smoothing;

  return this;
};

Material.prototype.setIsElastic = function(isElastic) {
  this.isElastic = isElastic;

  return this;
};

Material.prototype.setSpringK = function(springK) {
  this.springK = springK;

  return this;
};

Material.prototype.setYieldPoint = function(yieldPoint) {
  this.yieldPoint = yieldPoint;

  return this;
};

Material.prototype.setYieldRate = function(yieldRate) {
  this.yieldRate = yieldRate;

  return this;
};

var Node = function() {
	    this.mass = 0;
	    this.d = 0;
	    this.gx = 0;
	    this.gy = 0;
		// TODO: this currently limits the number of different materials that are available
	    this.cgx = [0, 0, 0, 0, 0, 0];
	    this.cgy = [0, 0, 0, 0, 0, 0];
	    this.velocity = Vector2.Zero.copy();
	    this.velocity2 = Vector2.Zero.copy();
	    this.acceleration = Vector2.Zero.copy();
	    
	    this.particleDensity = 0;
	};

var defaultNode = new Node();
	
	var Particle = function(x, y, u, v, material){
	    this.position = new Vector2(x, y);
	    this.prevPosition = new Vector2(x, y);
	    this.velocity = new Vector2(u, v);
	    // velocity gathered by the filter over the grid
	    this.gridVelocity = this.velocity.copy(); // or gradient x, y????

	    this.material = material;
	    
	    this.cellX = 0; // belongs to cell at x
	    this.cellY = 0; // belongs to cell at y

	    this.px = [0,0,0]; // deformation gradient?
	    this.py = [0,0,0];
	    this.gx = [0,0,0];
	    this.gy = [0,0,0];
	    
	    this.s = [0,0,0,0,0,0,0,0,0];
	    this.sx = [0,0,0,0,0,0,0,0,0];
	    this.sy = [0,0,0,0,0,0,0,0,0];
	    
	    this.node = [defaultNode, defaultNode, defaultNode, 
	                 defaultNode, defaultNode, defaultNode, 
	                 defaultNode, defaultNode, defaultNode];
	    
	    this.T00 = 0;
	    this.T01 = 0;
	    this.T11 = 0;
	};

var Obstacle = function(x, y, radius) {
		this.position = new Vector2(x, y);
		this.radius = radius;
	};

var Grid = function() {
  this.arr = [];
  this.activeCount = 0;
  this.gsizeY = 0;
  this.boundaries = new AABB();
  this.cellSize = Vector2.One.copy();
};

Grid.prototype.update = function(system) {
  this.recalculateBoundaries(system);
  this.clear();
  this.recalculateSizeY();
};

// TODO: reuse old grid
Grid.prototype.clear = function() {
  this.arr.length = 0;
  this.activeCount = 0;
};

Grid.prototype.iterate = function(fn, context) {
  var numberOfNodes = this.arr.length;
  for (var nIndex = 0; nIndex < numberOfNodes; nIndex++) {
    var n = this.arr[nIndex];
    if (n) {
      fn.call(context, n);
    }
  }
};

Grid.prototype.getAt = function(cellX, cellY) {
  return this.arr[cellX * this.gsizeY + cellY];
};

Grid.prototype.getOrCreateAt = function(cellX, cellY) {
  var cell = cellX * this.gsizeY + cellY;
  var node = this.arr[cell];

  if (node === undefined) {
    this.arr[cell] = node = new Node();
    this.activeCount++;
  }

  return node;
};

Grid.prototype.recalculateBoundaries = function(system) {
  // expand boundaries to include all particles
  this.boundaries.clear();
  system.particles.forEach(p => {
    this.boundaries.expandToInclude(p.position);
  });

  // expand boundaries a bit further
  this.boundaries.Min.x = Math.floor(this.boundaries.Min.x - 1);
  this.boundaries.Min.y = Math.floor(this.boundaries.Min.y - 1);
  this.boundaries.Max.x = Math.floor(this.boundaries.Max.x + 3);
  this.boundaries.Max.y = Math.floor(this.boundaries.Max.y + 3);
};

Grid.prototype.recalculateSizeY = function() {
  this.gsizeY = Math.floor(this.boundaries.Max.y - this.boundaries.Min.y);
};

var Integrator = function(grid) {
		this.grid = grid;
	};

	Integrator.prototype.updateStateAndGradientOf = function(particle) {
		var p = particle;
		// determine cell index for mesh
		p.cellX = Math.floor(p.position.x - this.grid.boundaries.Min.x - 0.5); // get cell x
		p.cellY = Math.floor(p.position.y - this.grid.boundaries.Min.y - 0.5); // get cell y

		var x = p.cellX - (p.position.x - this.grid.boundaries.Min.x);
		p.px[0] = (0.5 * x * x + 1.5 * x + 1.125);
		p.gx[0] = (x++ + 1.5);
		p.px[1] = (-x * x + 0.75);
		p.gx[1] = (-2.0 * (x++));
		p.px[2] = (0.5 * x * x - 1.5 * x + 1.125);
		p.gx[2] = (x - 1.5);

		var y = p.cellY - (p.position.y - this.grid.boundaries.Min.y);
		p.py[0] = (0.5 * y * y + 1.5 * y + 1.125);
		p.gy[0] = (y++ + 1.5);
		p.py[1] = (-y * y + 0.75);
		p.gy[1] = (-2.0 * (y++));
		p.py[2] = (0.5 * y * y - 1.5 * y + 1.125);
		p.gy[2] = (y - 1.5);

		// using quadratic interpolation
		// indices refer to corresponding adjacent cell
		
		// y  +-+-+-+
		//  2 |4|3|2|
		//    +-+-+-+
		//  1 |5|0|1|
		//    +-+-+-+
		//  0 |6|7|8|
		//    +-+-+-+
		//   /
		//  /  0 1 2 x
		
		// state variable
		p.s[0] = p.px[1] * p.py[1];
		p.s[1] = p.px[2] * p.py[1];
		p.s[2] = p.px[2] * p.py[2];
		p.s[3] = p.px[1] * p.py[2];
		p.s[4] = p.px[0] * p.py[2];
		p.s[5] = p.px[0] * p.py[1];
		p.s[6] = p.px[0] * p.py[0];
		p.s[7] = p.px[1] * p.py[0];
		p.s[8] = p.px[2] * p.py[0];
		
		// gradient in x axis
		p.sx[0] = p.gx[1] * p.py[1];
		p.sx[1] = p.gx[2] * p.py[1];
		p.sx[2] = p.gx[2] * p.py[2];
		p.sx[3] = p.gx[1] * p.py[2];
		p.sx[4] = p.gx[0] * p.py[2];
		p.sx[5] = p.gx[0] * p.py[1];
		p.sx[6] = p.gx[0] * p.py[0];
		p.sx[7] = p.gx[1] * p.py[0];
		p.sx[8] = p.gx[2] * p.py[0];

		// gradient in y axis
		p.sy[0] = p.px[1] * p.gy[1];
		p.sy[1] = p.px[2] * p.gy[1];
		p.sy[2] = p.px[2] * p.gy[2];
		p.sy[3] = p.px[1] * p.gy[2];
		p.sy[4] = p.px[0] * p.gy[2];
		p.sy[5] = p.px[0] * p.gy[1];
		p.sy[6] = p.px[0] * p.gy[0];
		p.sy[7] = p.px[1] * p.gy[0];
		p.sy[8] = p.px[2] * p.gy[0];
	};
	
	// cache the nodes a particle is adjacent to as a particles attribute
	Integrator.prototype.prepareParticle = function(particle) {
		particle.node[0] = this.grid.getOrCreateAt(particle.cellX+1, particle.cellY+1);
		particle.node[1] = this.grid.getOrCreateAt(particle.cellX+2, particle.cellY+1);
		particle.node[2] = this.grid.getOrCreateAt(particle.cellX+2, particle.cellY+2);
		particle.node[3] = this.grid.getOrCreateAt(particle.cellX+1, particle.cellY+2);
		particle.node[4] = this.grid.getOrCreateAt(particle.cellX,   particle.cellY+2);
		particle.node[5] = this.grid.getOrCreateAt(particle.cellX,   particle.cellY+1);
		particle.node[6] = this.grid.getOrCreateAt(particle.cellX,   particle.cellY  );
		particle.node[7] = this.grid.getOrCreateAt(particle.cellX+1, particle.cellY  );
		particle.node[8] = this.grid.getOrCreateAt(particle.cellX+2, particle.cellY  );
	};
	
	Integrator.prototype.integrate = function(particle, fn) {
		for(var i = 0; i < 9; i++)
			fn.call(undefined, particle, particle.node[i], particle.s[i], particle.sx[i], particle.sy[i]);
	};

var System = function() {
  this.wall = new AABB(new Vector2(-50, 2), new Vector2(50, 100));
  this.gravity = new Vector2(0, -0.05); // 0.004, 0.02
  this.materials = [];
  this.particles = [];
  this.springs = [];
  this.grid = new Grid();
  this.integrator = new Integrator(this.grid);

  this.useSurfaceTensionImplementation = true;
  this.drawGrid = false;

  this.doObstacles = false;
  this.obstacles = [];

  this.doSprings = false;
  this.drawSprings = false;
};

System.prototype.getNumberOfParticles = function() {
  return this.particles.length;
};

System.prototype.getNumberOfMaterials = function() {
  return this.materials.length;
};

System.prototype.createNewMaterial = function() {
  var newMaterial = new Material(this.materials.length);
  this.materials.push(newMaterial);

  return newMaterial;
};

System.prototype.addParticle = function(particle) {
  this.particles.push(particle);
};

/*
	 * UPDATE
	 */
System.prototype.update = function() {
  this.grid.update(this);

  if (this.useSurfaceTensionImplementation) {
    this.surfaceTensionImplementation();
  } else {
    this.simpleSimulation();
  }
};

/*
	 * surface tension implementation
	 */
System.prototype.surfaceTensionImplementation = function() {
  this.mapPropertiesToGrid();
  this.sumUpPerMaterialGradients();

  // Calculate pressure and add forces to grid
  this.particles.forEach((p, pIndex) => {
    var material = p.material;
    var dudx = 0,
      dudy = 0,
      dvdx = 0,
      dvdy = 0,
      sx = 0,
      sy = 0;

    this.integrator.integrate(p, function(particle, node, phi, gxpy, pxgy) {
      // Surface tension
      sx += phi * node.cgx[material.materialIndex];
      sy += phi * node.cgy[material.materialIndex];

      // Velocity gradient
      dudx += node.velocity.x * gxpy;
      dudy += node.velocity.x * pxgy;
      dvdx += node.velocity.y * gxpy;
      dvdy += node.velocity.y * pxgy;
    });

    // determine cell index for mesh
    var linearCellX = Math.floor(p.position.x - this.grid.boundaries.Min.x); // get cell x
    var linearCellY = Math.floor(p.position.y - this.grid.boundaries.Min.y); // get cell y

    // linear 2x2 kernel
    // y  +-+-+
    //  1 |2|4|
    //    +-+-+
    //  0 |1|3|
    //    +-+-+
    //   /
    //  /  0 1 x
    var n1 = this.grid.getOrCreateAt(linearCellX, linearCellY);
    var n2 = this.grid.getOrCreateAt(linearCellX, linearCellY + 1);
    var n3 = this.grid.getOrCreateAt(linearCellX + 1, linearCellY);
    var n4 = this.grid.getOrCreateAt(linearCellX + 1, linearCellY + 1);

    var density = this.uscip(
      n1.particleDensity,
      n1.gx,
      n1.gy,
      n2.particleDensity,
      n2.gx,
      n2.gy,
      n3.particleDensity,
      n3.gx,
      n3.gy,
      n4.particleDensity,
      n4.gx,
      n4.gy,
      p.position.x - this.grid.boundaries.Min.x - linearCellX,
      p.position.y - this.grid.boundaries.Min.y - linearCellY
    ); // r and s

    var restDensity = material.restDensity;
    var pressure = (material.stiffness * (density - restDensity)) / restDensity;
    if (pressure > 2.0) pressure = 2.0;

    // Update stress tensor
    var w1 = dudy - dvdx;
    var wT0 = 0.5 * w1 * (p.T01 + p.T01);
    var wT1 = 0.5 * w1 * (p.T00 - p.T11);
    var D00 = dudx;
    var D01 = 0.5 * (dudy + dvdx);
    var D11 = dvdy;
    var trace = 0.5 * (D00 + D11);
    p.T00 += 0.5 * (-wT0 + (D00 - trace) - material.meltRate * p.T00);
    p.T01 += 0.5 * (wT1 + D01 - material.meltRate * p.T01);
    p.T11 += 0.5 * (wT0 + (D11 - trace) - material.meltRate * p.T11);

    // Stress tensor fracture
    var norm = p.T00 * p.T00 + 2 * p.T01 * p.T01 + p.T11 * p.T11;

    if (norm > material.maxDeformation) {
      p.T00 = p.T01 = p.T11 = 0;
    }

    var T00 =
      material.particleMass *
      (material.elasticity * p.T00 +
        material.viscosity * D00 +
        pressure +
        trace * material.bulkViscosity);
    var T01 =
      material.particleMass *
      (material.elasticity * p.T01 + material.viscosity * D01);
    var T11 =
      material.particleMass *
      (material.elasticity * p.T11 +
        material.viscosity * D11 +
        pressure +
        trace * material.bulkViscosity);

    // Surface tension
    var lenSq = sx * sx + sy * sy;
    if (lenSq > 0) {
      var len = Math.sqrt(lenSq);
      var a = (material.particleMass * material.surfaceTension) / len;
      T00 -= a * (0.5 * lenSq - sx * sx);
      T01 -= a * (-sx * sy);
      T11 -= a * (0.5 * lenSq - sy * sy);
    }

    // Wall force
    var f = Vector2.Zero.copy();
    if (p.position.x < this.wall.Min.x) {
      f.x += this.wall.Min.x - p.position.x;
      p.velocity.x *= 0.1;
    }
    if (p.position.x > this.wall.Max.x) {
      f.x += this.wall.Max.x - p.position.x;
      p.velocity.x *= 0.1;
    }
    if (p.position.y < this.wall.Min.y) {
      f.y += this.wall.Min.y - p.position.y;
      p.velocity.y *= 0.1;
    }
    if (p.position.y > this.wall.Max.y) {
      f.y += this.wall.Max.y - p.position.y;
      p.velocity.y *= 0.1;
    }

    // test obstacle collision
    if (this.doObstacles) {
      // circular obstacles
      this.obstacles.forEach(obstacle => {
        var obstacleRadius = obstacle.radius;
        var obstacleRadiusSquared = obstacleRadius * obstacleRadius;
        var particleDistanceToMiddlePoint = obstacle.position.sub(p.position);
        var distanceSquared = particleDistanceToMiddlePoint.lengthSquared();
        if (distanceSquared < obstacleRadiusSquared) {
          var distance = Math.sqrt(distanceSquared);
          var dR = obstacleRadius - distance;
          f.subSelf(particleDistanceToMiddlePoint.mulFloat(dR / distance));
        }
      });
    }

    this.integrator.integrate(p, function(particle, node, phi, gxpy, pxgy) {
      node.acceleration.x += -(gxpy * T00 + pxgy * T01) + f.x * phi;
      node.acceleration.y += -(gxpy * T01 + pxgy * T11) + f.y * phi;
    });
  }, this);

  // Update acceleration of nodes
  this.grid.iterate(function(node) {
    if (node.mass > 0.0) {
      node.acceleration.divFloatSelf(node.mass);
    }
  }, this);

  this.particles.forEach((p, pIndex) => {
    var material = p.material;

    // Update particle velocities
    this.integrator.integrate(p, function(particle, node, phi, gxpy, pxgy) {
      particle.velocity.weightedAddSelf(node.acceleration, phi);
    });

    p.velocity.addSelf(this.gravity);
    p.velocity.mulFloatSelf(1 - material.damping);

    // Add particle velocities back to the grid
    this.integrator.integrate(p, function(particle, node, phi, gxpy, pxgy) {
      node.velocity2.weightedAddSelf(
        particle.velocity,
        material.particleMass * phi
      );
    });
  });

  // Update node velocities
  this.grid.iterate(function(node) {
    if (node.mass > 0.0) {
      node.velocity2.divFloatSelf(node.mass);
    }
  }, this);

  this.advanceParticles();
  this.springDisplacement();
  this.boundaryCorrection();
};

System.prototype.mapPropertiesToGrid = function() {
  this.particles.forEach((p, pIndex) => {
    var material = p.material;

    // Update grid cell index and kernel weights
    this.integrator.updateStateAndGradientOf(p);
    this.integrator.prepareParticle(p);

    // Add particle mass, velocity and density gradient to grid
    this.integrator.integrate(p, function(particle, node, phi, gxpy, pxgy) {
      node.mass += phi * material.particleMass;
      node.particleDensity += phi;
      node.velocity.weightedAddSelf(
        particle.velocity,
        material.particleMass * phi
      );
      node.cgx[material.materialIndex] += gxpy;
      node.cgy[material.materialIndex] += pxgy;
    });
  });
};

System.prototype.sumUpPerMaterialGradients = function() {
  var numMaterials = this.getNumberOfMaterials();
  this.grid.iterate(function(node) {
    if (node.mass > 0.0) {
      node.acceleration.clear();
      node.gx = 0;
      node.gy = 0;
      node.velocity.divFloatSelf(node.mass);
      // sum up gradients of all materials
      for (var i = 0; i < numMaterials; i++) {
        node.gx += node.cgx[i];
        node.gy += node.cgy[i];
      }
      for (var i = 0; i < numMaterials; i++) {
        node.cgx[i] -= node.gx - node.cgx[i];
        node.cgy[i] -= node.gy - node.cgy[i];
      }
    }
  }, this);
};

System.prototype.uscip = function(
  p00,
  x00,
  y00,
  p01,
  x01,
  y01,
  p10,
  x10,
  y10,
  p11,
  x11,
  y11,
  u,
  v
) {
  var dx = x00 - x01;
  var dy = y00 - y10;
  var a = p01 - p00;
  var b = p11 - p10 - a;
  var c = p10 - p00;
  var d = y11 - y01;
  return (
    ((((d - 2 * b - dy) * u - 2 * a + y00 + y01) * v +
      ((3 * b + 2 * dy - d) * u + 3 * a - 2 * y00 - y01)) *
      v +
      ((((2 * c - x00 - x10) * u + (3 * b + 2 * dx + x10 - x11)) * u -
        b -
        dy -
        dx) *
        u +
        y00)) *
      v +
    (((x11 - 2 * (p11 - p01 + c) + x10 + x00 + x01) * u +
      (3 * c - 2 * x00 - x10)) *
      u +
      x00) *
      u +
    p00
  );
};

System.prototype.advanceParticles = function() {
  this.particles.forEach((p, pIndex) => {
    var material = p.material;

    var gVelocity = Vector2.Zero.copy();
    var dudx = 0,
      dudy = 0,
      dvdx = 0,
      dvdy = 0;

    this.integrator.integrate(p, function(particle, node, phi, gxpy, pxgy) {
      gVelocity.weightedAddSelf(node.velocity2, phi);

      // Velocity gradient
      dudx += node.velocity2.x * gxpy;
      dudy += node.velocity2.x * pxgy;
      dvdx += node.velocity2.y * gxpy;
      dvdy += node.velocity2.y * pxgy;
    });

    // Update stress tensor
    var w1 = dudy - dvdx;
    var wT0 = 0.5 * w1 * (p.T01 + p.T01);
    var wT1 = 0.5 * w1 * (p.T00 - p.T11);
    var D00 = dudx;
    var D01 = 0.5 * (dudy + dvdx);
    var D11 = dvdy;
    var trace = 0.5 * (D00 + D11);
    p.T00 += 0.5 * (-wT0 + (D00 - trace) - material.meltRate * p.T00);
    p.T01 += 0.5 * (wT1 + D01 - material.meltRate * p.T01);
    p.T11 += 0.5 * (wT0 + (D11 - trace) - material.meltRate * p.T11);

    // Stress tensor fracture
    var norm = p.T00 * p.T00 + 2 * p.T01 * p.T01 + p.T11 * p.T11;

    if (norm > material.maxDeformation) {
      p.T00 = p.T01 = p.T11 = 0;
    }

    // advance particle
    p.prevPosition.set(p.position);
    p.position.addSelf(gVelocity);
    p.gridVelocity.set(gVelocity);
    p.velocity.weightedAddSelf(gVelocity.sub(p.velocity), material.smoothing);
  });
};

System.prototype.springDisplacement = function() {
  if (this.doSprings) {
    this.springs.forEach((s, sIndex) => {
      s.update();
      s.solve();
    });
  }
};

// hard boundary correction
System.prototype.boundaryCorrection = function() {
  this.particles.forEach((p, pIndex) => {
    if (p.position.x < this.wall.Min.x - 4)
      p.position.x = this.wall.Min.x - 4 + 0.01 * Math.random();
    else if (p.position.x > this.wall.Max.x + 4)
      p.position.x = this.wall.Max.x + 4 - 0.01 * Math.random();
    if (p.position.y < this.wall.Min.y - 4)
      p.position.y = this.wall.Min.y - 4 + 0.01 * Math.random();
    else if (p.position.y > this.wall.Max.y + 4)
      p.position.y = this.wall.Max.y + 4 - 0.01 * Math.random();
  });
};

var Spring = function(particle1, particle2, restLength) {
		this.particle1 = particle1;
		this.particle2 = particle2;
		this.restLength = restLength;
		this.currentDistance = 0;
	};
	
	Spring.prototype.update = function() {
		this.currentDistance = this.particle1.position.distance(this.particle2.position);
	};
	
	Spring.prototype.solve = function() {
		var p2position = this.particle2.position;
		var p1position = this.particle1.position;
		var rij = p2position.sub(p1position);
		rij.mulFloatSelf(1 / this.currentDistance);
		var D = this.particle1.material.springK * (this.restLength - this.currentDistance);
		rij.mulFloatSelf(D * 0.5);
		p1position.subSelf(rij);
		p2position.addSelf(rij);
	};
	
	Spring.prototype.contains = function(particle) {
		return this.particle1 === particle || this.particle2 === particle;
	};

var Group = function(system, minX, minY, maxX, maxY, u, v, material) {
		this.material = material;
		
		var map = [];
		for (var i = minX; i < maxX; i++) {
			map[map.length] = [];
	        for (var j = minY; j < maxY; j++) {
	        	var p = new Particle(i, j, u, v, material);
	        	system.addParticle(p);
	        	if(material.isElastic) {
	        		map[map.length-1].push(p);
	        	}
	        }
		}
    	if(material.isElastic) {
    		for(var i = 0; i < map.length; i++) {
    			for(var j = 1; j < map[0].length; j++) {
        			system.springs.push(new Spring(map[i][j-1], map[i][j], map[i][j-1].position.distance(map[i][j].position)));
    			}
    		}
    		for(var j = 0; j < map[0].length; j++) {
    			for(var i = 1; i < map.length; i++) {
        			system.springs.push(new Spring(map[i-1][j], map[i][j], map[i-1][j].position.distance(map[i][j].position)));
    			}
    		}
    		for(var i = 1; i < map.length; i++) {
    			for(var j = 1; j < map[0].length; j++) {
        			system.springs.push(new Spring(map[i-1][j-1], map[i][j], map[i-1][j-1].position.distance(map[i][j].position)));
    			}
    		}
    		for(var i = 0; i < map.length - 1; i++) {
    			for(var j = 1; j < map[0].length; j++) {
        			system.springs.push(new Spring(map[i+1][j-1], map[i][j], map[i+1][j-1].position.distance(map[i][j].position)));
    			}
    		}
    	}
	};

/*
	 * Early simple implementation of the material point method
	 */
System.prototype.simpleSimulation = function() {
  this.__calculateParticleKernels();
  this.__sumParticleDensityFromGridAndAddPressureansElasticForcesToGrid();
  this.__divideGridAccelerationByMass();
  this.__accelerateParticlesAndInterpolateVelocityBackToGrid();
  this.__divideGridVelocityByMass();
  this.__advanceParticles();
};

System.prototype.__calculateParticleKernels = function() {
  // calculate particle kernels, and add density and density gradients to the grid
  this.particles.forEach((p, pIndex) => {
    this.integrator.updateStateAndGradientOf(p);
    this.integrator.prepareParticle(p);

    this.integrator.integrate(p, function(particle, node, phi, gxpy, pxgy) {
      node.mass += phi;
    });
  });
};

System.prototype.__sumParticleDensityFromGridAndAddPressureansElasticForcesToGrid = function() {
  // Sum particle density from grid, and add pressure and elastic forces to grid
  this.particles.forEach((p, pIndex) => {
    var density = 0;
    this.integrator.integrate(p, function(particle, node, phi, gxpy, pxgy) {
      density += phi * node.mass;
    });

    var restDensity = p.material.restDensity;
    var pressure = (density - restDensity) / restDensity;
    if (pressure > 4.0) pressure = 4.0;

    var f = Vector2.Zero.copy();
    if (p.position.x < this.wall.Min.x) {
      f.x += this.wall.Min.x - p.position.x;
      p.velocity.x *= 0.1;
    }
    if (p.position.x > this.wall.Max.x) {
      f.x += this.wall.Max.x - p.position.x;
      p.velocity.x *= 0.1;
    }
    if (p.position.y < this.wall.Min.y) {
      f.y += this.wall.Min.y - p.position.y;
      p.velocity.y *= 0.1;
    }
    if (p.position.y > this.wall.Max.y) {
      f.y += this.wall.Max.y - p.position.y;
      p.velocity.y *= 0.1;
    }

    this.integrator.integrate(p, function(particle, node, phi, gxpy, pxgy) {
      node.acceleration.x += -(gxpy * pressure) + f.x * phi;
      node.acceleration.y += -(pxgy * pressure) + f.y * phi;
    });
  });
};

// divide grid acceleration by mass
System.prototype.__divideGridAccelerationByMass = function() {
  this.grid.iterate(function(n) {
    if (n.mass > 0.0) {
      n.acceleration.divFloatSelf(n.mass);
      n.acceleration.addSelf(this.gravity);
    }
  }, this);
};

// accelerate particles and interpolate velocity back to grid
System.prototype.__accelerateParticlesAndInterpolateVelocityBackToGrid = function() {
  this.particles.forEach((p, pIndex) => {
    this.integrator.integrate(p, function(particle, node, phi, gxpy, pxgy) {
      particle.velocity.weightedAddSelf(node.acceleration, phi);
    });
    this.integrator.integrate(p, function(particle, node, phi, gxpy, pxgy) {
      node.velocity.weightedAddSelf(particle.velocity, phi);
    });
  });
};

// divide grid velocity by mass
System.prototype.__divideGridVelocityByMass = function() {
  this.grid.iterate(function(n) {
    if (n.mass > 0.0) {
      n.velocity.divFloatSelf(n.mass);
    }
  }, this);
};

System.prototype.__advanceParticles = function() {
  // advance particles
  this.particles.forEach((p, pIndex) => {
    p.gridVelocity.clear();
    this.integrator.integrate(p, function(particle, node, phi, gxpy, pxgy) {
      particle.gridVelocity.weightedAddSelf(node.velocity, phi);
    });
    p.position.addSelf(p.gridVelocity);
    p.velocity.weightedAddSelf(
      p.gridVelocity.sub(p.velocity),
      p.material.smoothing
    );
  });
};

var Floom = function() {};

Floom.Material = Material;
Floom.Particle = Particle;
Floom.Group = Group;
Floom.Node = Node;
Floom.Spring = Spring;
Floom.Obstacle = Obstacle;
Floom.System = System;

// export { default as Floom } from "./floom/floom.js";

export default Floom;
export { AABB, Vector2 };
