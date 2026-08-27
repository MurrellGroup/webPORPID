//#region src/binary.ts
const encoder$1 = new TextEncoder();
const decoder$1 = new TextDecoder();
var BinaryWriter = class {
	chunks = [];
	length = 0;
	push(bytes) {
		this.chunks.push(bytes);
		this.length += bytes.byteLength;
	}
	magic(value) {
		this.push(encoder$1.encode(value));
	}
	u8(value) {
		this.push(Uint8Array.of(value & 255));
	}
	u16(value) {
		const b = new Uint8Array(2);
		new DataView(b.buffer).setUint16(0, value, true);
		this.push(b);
	}
	i32(value) {
		const b = new Uint8Array(4);
		new DataView(b.buffer).setInt32(0, value, true);
		this.push(b);
	}
	u32(value) {
		const b = new Uint8Array(4);
		new DataView(b.buffer).setUint32(0, value, true);
		this.push(b);
	}
	u64(value) {
		const b = new Uint8Array(8);
		new DataView(b.buffer).setBigUint64(0, value, true);
		this.push(b);
	}
	f32(value) {
		const b = new Uint8Array(4);
		new DataView(b.buffer).setFloat32(0, value, true);
		this.push(b);
	}
	f64(value) {
		const b = new Uint8Array(8);
		new DataView(b.buffer).setFloat64(0, value, true);
		this.push(b);
	}
	bytes(value) {
		this.push(value);
	}
	string(value) {
		const bytes = encoder$1.encode(value);
		this.u32(bytes.byteLength);
		this.push(bytes);
	}
	finish() {
		const output = new Uint8Array(this.length);
		let offset = 0;
		for (const chunk of this.chunks) {
			output.set(chunk, offset);
			offset += chunk.byteLength;
		}
		return output;
	}
};
var BinaryReader = class {
	bytesValue;
	view;
	offset = 0;
	constructor(bytesValue) {
		this.bytesValue = bytesValue;
		this.view = new DataView(bytesValue.buffer, bytesValue.byteOffset, bytesValue.byteLength);
	}
	require(length) {
		if (this.offset + length > this.bytesValue.byteLength) throw new Error("Truncated webPORPID binary payload.");
	}
	magic(value) {
		const expected = encoder$1.encode(value);
		this.require(expected.length);
		for (let index = 0; index < expected.length; index++) if (this.bytesValue[this.offset + index] !== expected[index]) throw new Error(`Expected ${value} payload.`);
		this.offset += expected.length;
	}
	u8() {
		this.require(1);
		return this.view.getUint8(this.offset++);
	}
	u16() {
		this.require(2);
		const v = this.view.getUint16(this.offset, true);
		this.offset += 2;
		return v;
	}
	i32() {
		this.require(4);
		const v = this.view.getInt32(this.offset, true);
		this.offset += 4;
		return v;
	}
	u32() {
		this.require(4);
		const v = this.view.getUint32(this.offset, true);
		this.offset += 4;
		return v;
	}
	u64() {
		this.require(8);
		const v = this.view.getBigUint64(this.offset, true);
		this.offset += 8;
		return v;
	}
	f32() {
		this.require(4);
		const v = this.view.getFloat32(this.offset, true);
		this.offset += 4;
		return v;
	}
	f64() {
		this.require(8);
		const v = this.view.getFloat64(this.offset, true);
		this.offset += 8;
		return v;
	}
	bytes(length) {
		this.require(length);
		const v = this.bytesValue.subarray(this.offset, this.offset + length);
		this.offset += length;
		return v;
	}
	string() {
		const length = this.u32();
		return decoder$1.decode(this.bytes(length));
	}
	get done() {
		return this.offset === this.bytesValue.byteLength;
	}
	get remaining() {
		return this.bytesValue.byteLength - this.offset;
	}
};
var Iovec = class Iovec {
	static read_bytes(view, ptr) {
		const iovec = new Iovec();
		iovec.buf = view.getUint32(ptr, true);
		iovec.buf_len = view.getUint32(ptr + 4, true);
		return iovec;
	}
	static read_bytes_array(view, ptr, len) {
		const iovecs = [];
		for (let i = 0; i < len; i++) iovecs.push(Iovec.read_bytes(view, ptr + 8 * i));
		return iovecs;
	}
};
var Ciovec = class Ciovec {
	static read_bytes(view, ptr) {
		const iovec = new Ciovec();
		iovec.buf = view.getUint32(ptr, true);
		iovec.buf_len = view.getUint32(ptr + 4, true);
		return iovec;
	}
	static read_bytes_array(view, ptr, len) {
		const iovecs = [];
		for (let i = 0; i < len; i++) iovecs.push(Ciovec.read_bytes(view, ptr + 8 * i));
		return iovecs;
	}
};
var Subscription = class Subscription {
	static read_bytes(view, ptr) {
		return new Subscription(view.getBigUint64(ptr, true), view.getUint8(ptr + 8), view.getUint32(ptr + 16, true), view.getBigUint64(ptr + 24, true), view.getUint16(ptr + 36, true));
	}
	constructor(userdata, eventtype, clockid, timeout, flags) {
		this.userdata = userdata;
		this.eventtype = eventtype;
		this.clockid = clockid;
		this.timeout = timeout;
		this.flags = flags;
	}
};
var Event = class {
	write_bytes(view, ptr) {
		view.setBigUint64(ptr, this.userdata, true);
		view.setUint16(ptr + 8, this.error, true);
		view.setUint8(ptr + 10, this.eventtype);
	}
	constructor(userdata, error, eventtype) {
		this.userdata = userdata;
		this.error = error;
		this.eventtype = eventtype;
	}
};
//#endregion
//#region node_modules/@bjorn3/browser_wasi_shim/dist/debug.js
let Debug = class Debug {
	enable(enabled) {
		this.log = createLogger(enabled === void 0 ? true : enabled, this.prefix);
	}
	get enabled() {
		return this.isEnabled;
	}
	constructor(isEnabled) {
		this.isEnabled = isEnabled;
		this.prefix = "wasi:";
		this.enable(isEnabled);
	}
};
function createLogger(enabled, prefix) {
	if (enabled) return console.log.bind(console, "%c%s", "color: #265BA0", prefix);
	else return () => {};
}
const debug = new Debug(false);
//#endregion
//#region node_modules/@bjorn3/browser_wasi_shim/dist/wasi.js
var WASIProcExit = class extends Error {
	constructor(code) {
		super("exit with exit code " + code);
		this.code = code;
	}
};
let WASI = class WASI {
	start(instance) {
		this.inst = instance;
		try {
			instance.exports._start();
			return 0;
		} catch (e) {
			if (e instanceof WASIProcExit) return e.code;
			else throw e;
		}
	}
	initialize(instance) {
		this.inst = instance;
		if (instance.exports._initialize) instance.exports._initialize();
	}
	constructor(args, env, fds, options = {}) {
		this.args = [];
		this.env = [];
		this.fds = [];
		debug.enable(options.debug);
		this.args = args;
		this.env = env;
		this.fds = fds;
		const self = this;
		this.wasiImport = {
			args_sizes_get(argc, argv_buf_size) {
				const buffer = new DataView(self.inst.exports.memory.buffer);
				buffer.setUint32(argc, self.args.length, true);
				let buf_size = 0;
				for (const arg of self.args) buf_size += arg.length + 1;
				buffer.setUint32(argv_buf_size, buf_size, true);
				debug.log(buffer.getUint32(argc, true), buffer.getUint32(argv_buf_size, true));
				return 0;
			},
			args_get(argv, argv_buf) {
				const buffer = new DataView(self.inst.exports.memory.buffer);
				const buffer8 = new Uint8Array(self.inst.exports.memory.buffer);
				const orig_argv_buf = argv_buf;
				for (let i = 0; i < self.args.length; i++) {
					buffer.setUint32(argv, argv_buf, true);
					argv += 4;
					const arg = new TextEncoder().encode(self.args[i]);
					buffer8.set(arg, argv_buf);
					buffer.setUint8(argv_buf + arg.length, 0);
					argv_buf += arg.length + 1;
				}
				if (debug.enabled) debug.log(new TextDecoder("utf-8").decode(buffer8.slice(orig_argv_buf, argv_buf)));
				return 0;
			},
			environ_sizes_get(environ_count, environ_size) {
				const buffer = new DataView(self.inst.exports.memory.buffer);
				buffer.setUint32(environ_count, self.env.length, true);
				let buf_size = 0;
				for (const environ of self.env) buf_size += new TextEncoder().encode(environ).length + 1;
				buffer.setUint32(environ_size, buf_size, true);
				debug.log(buffer.getUint32(environ_count, true), buffer.getUint32(environ_size, true));
				return 0;
			},
			environ_get(environ, environ_buf) {
				const buffer = new DataView(self.inst.exports.memory.buffer);
				const buffer8 = new Uint8Array(self.inst.exports.memory.buffer);
				const orig_environ_buf = environ_buf;
				for (let i = 0; i < self.env.length; i++) {
					buffer.setUint32(environ, environ_buf, true);
					environ += 4;
					const e = new TextEncoder().encode(self.env[i]);
					buffer8.set(e, environ_buf);
					buffer.setUint8(environ_buf + e.length, 0);
					environ_buf += e.length + 1;
				}
				if (debug.enabled) debug.log(new TextDecoder("utf-8").decode(buffer8.slice(orig_environ_buf, environ_buf)));
				return 0;
			},
			clock_res_get(id, res_ptr) {
				let resolutionValue;
				switch (id) {
					case 1:
						resolutionValue = 5000n;
						break;
					case 0:
						resolutionValue = 1000000n;
						break;
					default: return 52;
				}
				new DataView(self.inst.exports.memory.buffer).setBigUint64(res_ptr, resolutionValue, true);
				return 0;
			},
			clock_time_get(id, precision, time) {
				const buffer = new DataView(self.inst.exports.memory.buffer);
				if (id === 0) buffer.setBigUint64(time, BigInt((/* @__PURE__ */ new Date()).getTime()) * 1000000n, true);
				else if (id == 1) {
					let monotonic_time;
					try {
						monotonic_time = BigInt(Math.round(performance.now() * 1e6));
					} catch (e) {
						monotonic_time = 0n;
					}
					buffer.setBigUint64(time, monotonic_time, true);
				} else buffer.setBigUint64(time, 0n, true);
				return 0;
			},
			fd_advise(fd, offset, len, advice) {
				if (self.fds[fd] != void 0) return 0;
				else return 8;
			},
			fd_allocate(fd, offset, len) {
				if (self.fds[fd] != void 0) return self.fds[fd].fd_allocate(offset, len);
				else return 8;
			},
			fd_close(fd) {
				if (self.fds[fd] != void 0) {
					const ret = self.fds[fd].fd_close();
					self.fds[fd] = void 0;
					return ret;
				} else return 8;
			},
			fd_datasync(fd) {
				if (self.fds[fd] != void 0) return self.fds[fd].fd_sync();
				else return 8;
			},
			fd_fdstat_get(fd, fdstat_ptr) {
				if (self.fds[fd] != void 0) {
					const { ret, fdstat } = self.fds[fd].fd_fdstat_get();
					if (fdstat != null) fdstat.write_bytes(new DataView(self.inst.exports.memory.buffer), fdstat_ptr);
					return ret;
				} else return 8;
			},
			fd_fdstat_set_flags(fd, flags) {
				if (self.fds[fd] != void 0) return self.fds[fd].fd_fdstat_set_flags(flags);
				else return 8;
			},
			fd_fdstat_set_rights(fd, fs_rights_base, fs_rights_inheriting) {
				if (self.fds[fd] != void 0) return self.fds[fd].fd_fdstat_set_rights(fs_rights_base, fs_rights_inheriting);
				else return 8;
			},
			fd_filestat_get(fd, filestat_ptr) {
				if (self.fds[fd] != void 0) {
					const { ret, filestat } = self.fds[fd].fd_filestat_get();
					if (filestat != null) filestat.write_bytes(new DataView(self.inst.exports.memory.buffer), filestat_ptr);
					return ret;
				} else return 8;
			},
			fd_filestat_set_size(fd, size) {
				if (self.fds[fd] != void 0) return self.fds[fd].fd_filestat_set_size(size);
				else return 8;
			},
			fd_filestat_set_times(fd, atim, mtim, fst_flags) {
				if (self.fds[fd] != void 0) return self.fds[fd].fd_filestat_set_times(atim, mtim, fst_flags);
				else return 8;
			},
			fd_pread(fd, iovs_ptr, iovs_len, offset, nread_ptr) {
				const buffer = new DataView(self.inst.exports.memory.buffer);
				const buffer8 = new Uint8Array(self.inst.exports.memory.buffer);
				if (self.fds[fd] != void 0) {
					const iovecs = Iovec.read_bytes_array(buffer, iovs_ptr, iovs_len);
					let nread = 0;
					for (const iovec of iovecs) {
						const { ret, data } = self.fds[fd].fd_pread(iovec.buf_len, offset);
						if (ret != 0) {
							buffer.setUint32(nread_ptr, nread, true);
							return ret;
						}
						buffer8.set(data, iovec.buf);
						nread += data.length;
						offset += BigInt(data.length);
						if (data.length != iovec.buf_len) break;
					}
					buffer.setUint32(nread_ptr, nread, true);
					return 0;
				} else return 8;
			},
			fd_prestat_get(fd, buf_ptr) {
				const buffer = new DataView(self.inst.exports.memory.buffer);
				if (self.fds[fd] != void 0) {
					const { ret, prestat } = self.fds[fd].fd_prestat_get();
					if (prestat != null) prestat.write_bytes(buffer, buf_ptr);
					return ret;
				} else return 8;
			},
			fd_prestat_dir_name(fd, path_ptr, path_len) {
				if (self.fds[fd] != void 0) {
					const { ret, prestat } = self.fds[fd].fd_prestat_get();
					if (prestat == null) return ret;
					const prestat_dir_name = prestat.inner.pr_name;
					new Uint8Array(self.inst.exports.memory.buffer).set(prestat_dir_name.slice(0, path_len), path_ptr);
					return prestat_dir_name.byteLength > path_len ? 37 : 0;
				} else return 8;
			},
			fd_pwrite(fd, iovs_ptr, iovs_len, offset, nwritten_ptr) {
				const buffer = new DataView(self.inst.exports.memory.buffer);
				const buffer8 = new Uint8Array(self.inst.exports.memory.buffer);
				if (self.fds[fd] != void 0) {
					const iovecs = Ciovec.read_bytes_array(buffer, iovs_ptr, iovs_len);
					let nwritten = 0;
					for (const iovec of iovecs) {
						const data = buffer8.slice(iovec.buf, iovec.buf + iovec.buf_len);
						const { ret, nwritten: nwritten_part } = self.fds[fd].fd_pwrite(data, offset);
						if (ret != 0) {
							buffer.setUint32(nwritten_ptr, nwritten, true);
							return ret;
						}
						nwritten += nwritten_part;
						offset += BigInt(nwritten_part);
						if (nwritten_part != data.byteLength) break;
					}
					buffer.setUint32(nwritten_ptr, nwritten, true);
					return 0;
				} else return 8;
			},
			fd_read(fd, iovs_ptr, iovs_len, nread_ptr) {
				const buffer = new DataView(self.inst.exports.memory.buffer);
				const buffer8 = new Uint8Array(self.inst.exports.memory.buffer);
				if (self.fds[fd] != void 0) {
					const iovecs = Iovec.read_bytes_array(buffer, iovs_ptr, iovs_len);
					let nread = 0;
					for (const iovec of iovecs) {
						const { ret, data } = self.fds[fd].fd_read(iovec.buf_len);
						if (ret != 0) {
							buffer.setUint32(nread_ptr, nread, true);
							return ret;
						}
						buffer8.set(data, iovec.buf);
						nread += data.length;
						if (data.length != iovec.buf_len) break;
					}
					buffer.setUint32(nread_ptr, nread, true);
					return 0;
				} else return 8;
			},
			fd_readdir(fd, buf, buf_len, cookie, bufused_ptr) {
				const buffer = new DataView(self.inst.exports.memory.buffer);
				const buffer8 = new Uint8Array(self.inst.exports.memory.buffer);
				if (self.fds[fd] != void 0) {
					let bufused = 0;
					while (true) {
						const { ret, dirent } = self.fds[fd].fd_readdir_single(cookie);
						if (ret != 0) {
							buffer.setUint32(bufused_ptr, bufused, true);
							return ret;
						}
						if (dirent == null) break;
						if (buf_len - bufused < dirent.head_length()) {
							bufused = buf_len;
							break;
						}
						const head_bytes = new ArrayBuffer(dirent.head_length());
						dirent.write_head_bytes(new DataView(head_bytes), 0);
						buffer8.set(new Uint8Array(head_bytes).slice(0, Math.min(head_bytes.byteLength, buf_len - bufused)), buf);
						buf += dirent.head_length();
						bufused += dirent.head_length();
						if (buf_len - bufused < dirent.name_length()) {
							bufused = buf_len;
							break;
						}
						dirent.write_name_bytes(buffer8, buf, buf_len - bufused);
						buf += dirent.name_length();
						bufused += dirent.name_length();
						cookie = dirent.d_next;
					}
					buffer.setUint32(bufused_ptr, bufused, true);
					return 0;
				} else return 8;
			},
			fd_renumber(fd, to) {
				if (self.fds[fd] != void 0 && self.fds[to] != void 0) {
					const ret = self.fds[to].fd_close();
					if (ret != 0) return ret;
					self.fds[to] = self.fds[fd];
					self.fds[fd] = void 0;
					return 0;
				} else return 8;
			},
			fd_seek(fd, offset, whence, offset_out_ptr) {
				const buffer = new DataView(self.inst.exports.memory.buffer);
				if (self.fds[fd] != void 0) {
					const { ret, offset: offset_out } = self.fds[fd].fd_seek(offset, whence);
					buffer.setBigInt64(offset_out_ptr, offset_out, true);
					return ret;
				} else return 8;
			},
			fd_sync(fd) {
				if (self.fds[fd] != void 0) return self.fds[fd].fd_sync();
				else return 8;
			},
			fd_tell(fd, offset_ptr) {
				const buffer = new DataView(self.inst.exports.memory.buffer);
				if (self.fds[fd] != void 0) {
					const { ret, offset } = self.fds[fd].fd_tell();
					buffer.setBigUint64(offset_ptr, offset, true);
					return ret;
				} else return 8;
			},
			fd_write(fd, iovs_ptr, iovs_len, nwritten_ptr) {
				const buffer = new DataView(self.inst.exports.memory.buffer);
				const buffer8 = new Uint8Array(self.inst.exports.memory.buffer);
				if (self.fds[fd] != void 0) {
					const iovecs = Ciovec.read_bytes_array(buffer, iovs_ptr, iovs_len);
					let nwritten = 0;
					for (const iovec of iovecs) {
						const data = buffer8.slice(iovec.buf, iovec.buf + iovec.buf_len);
						const { ret, nwritten: nwritten_part } = self.fds[fd].fd_write(data);
						if (ret != 0) {
							buffer.setUint32(nwritten_ptr, nwritten, true);
							return ret;
						}
						nwritten += nwritten_part;
						if (nwritten_part != data.byteLength) break;
					}
					buffer.setUint32(nwritten_ptr, nwritten, true);
					return 0;
				} else return 8;
			},
			path_create_directory(fd, path_ptr, path_len) {
				const buffer8 = new Uint8Array(self.inst.exports.memory.buffer);
				if (self.fds[fd] != void 0) {
					const path = new TextDecoder("utf-8").decode(buffer8.slice(path_ptr, path_ptr + path_len));
					return self.fds[fd].path_create_directory(path);
				} else return 8;
			},
			path_filestat_get(fd, flags, path_ptr, path_len, filestat_ptr) {
				const buffer = new DataView(self.inst.exports.memory.buffer);
				const buffer8 = new Uint8Array(self.inst.exports.memory.buffer);
				if (self.fds[fd] != void 0) {
					const path = new TextDecoder("utf-8").decode(buffer8.slice(path_ptr, path_ptr + path_len));
					const { ret, filestat } = self.fds[fd].path_filestat_get(flags, path);
					if (filestat != null) filestat.write_bytes(buffer, filestat_ptr);
					return ret;
				} else return 8;
			},
			path_filestat_set_times(fd, flags, path_ptr, path_len, atim, mtim, fst_flags) {
				const buffer8 = new Uint8Array(self.inst.exports.memory.buffer);
				if (self.fds[fd] != void 0) {
					const path = new TextDecoder("utf-8").decode(buffer8.slice(path_ptr, path_ptr + path_len));
					return self.fds[fd].path_filestat_set_times(flags, path, atim, mtim, fst_flags);
				} else return 8;
			},
			path_link(old_fd, old_flags, old_path_ptr, old_path_len, new_fd, new_path_ptr, new_path_len) {
				const buffer8 = new Uint8Array(self.inst.exports.memory.buffer);
				if (self.fds[old_fd] != void 0 && self.fds[new_fd] != void 0) {
					const old_path = new TextDecoder("utf-8").decode(buffer8.slice(old_path_ptr, old_path_ptr + old_path_len));
					const new_path = new TextDecoder("utf-8").decode(buffer8.slice(new_path_ptr, new_path_ptr + new_path_len));
					const { ret, inode_obj } = self.fds[old_fd].path_lookup(old_path, old_flags);
					if (inode_obj == null) return ret;
					return self.fds[new_fd].path_link(new_path, inode_obj, false);
				} else return 8;
			},
			path_open(fd, dirflags, path_ptr, path_len, oflags, fs_rights_base, fs_rights_inheriting, fd_flags, opened_fd_ptr) {
				const buffer = new DataView(self.inst.exports.memory.buffer);
				const buffer8 = new Uint8Array(self.inst.exports.memory.buffer);
				if (self.fds[fd] != void 0) {
					const path = new TextDecoder("utf-8").decode(buffer8.slice(path_ptr, path_ptr + path_len));
					debug.log(path);
					const { ret, fd_obj } = self.fds[fd].path_open(dirflags, path, oflags, fs_rights_base, fs_rights_inheriting, fd_flags);
					if (ret != 0) return ret;
					self.fds.push(fd_obj);
					const opened_fd = self.fds.length - 1;
					buffer.setUint32(opened_fd_ptr, opened_fd, true);
					return 0;
				} else return 8;
			},
			path_readlink(fd, path_ptr, path_len, buf_ptr, buf_len, nread_ptr) {
				const buffer = new DataView(self.inst.exports.memory.buffer);
				const buffer8 = new Uint8Array(self.inst.exports.memory.buffer);
				if (self.fds[fd] != void 0) {
					const path = new TextDecoder("utf-8").decode(buffer8.slice(path_ptr, path_ptr + path_len));
					debug.log(path);
					const { ret, data } = self.fds[fd].path_readlink(path);
					if (data != null) {
						const data_buf = new TextEncoder().encode(data);
						if (data_buf.length > buf_len) {
							buffer.setUint32(nread_ptr, 0, true);
							return 8;
						}
						buffer8.set(data_buf, buf_ptr);
						buffer.setUint32(nread_ptr, data_buf.length, true);
					}
					return ret;
				} else return 8;
			},
			path_remove_directory(fd, path_ptr, path_len) {
				const buffer8 = new Uint8Array(self.inst.exports.memory.buffer);
				if (self.fds[fd] != void 0) {
					const path = new TextDecoder("utf-8").decode(buffer8.slice(path_ptr, path_ptr + path_len));
					return self.fds[fd].path_remove_directory(path);
				} else return 8;
			},
			path_rename(fd, old_path_ptr, old_path_len, new_fd, new_path_ptr, new_path_len) {
				const buffer8 = new Uint8Array(self.inst.exports.memory.buffer);
				if (self.fds[fd] != void 0 && self.fds[new_fd] != void 0) {
					const old_path = new TextDecoder("utf-8").decode(buffer8.slice(old_path_ptr, old_path_ptr + old_path_len));
					const new_path = new TextDecoder("utf-8").decode(buffer8.slice(new_path_ptr, new_path_ptr + new_path_len));
					let { ret, inode_obj } = self.fds[fd].path_unlink(old_path);
					if (inode_obj == null) return ret;
					ret = self.fds[new_fd].path_link(new_path, inode_obj, true);
					if (ret != 0) {
						if (self.fds[fd].path_link(old_path, inode_obj, true) != 0) throw "path_link should always return success when relinking an inode back to the original place";
					}
					return ret;
				} else return 8;
			},
			path_symlink(old_path_ptr, old_path_len, fd, new_path_ptr, new_path_len) {
				const buffer8 = new Uint8Array(self.inst.exports.memory.buffer);
				if (self.fds[fd] != void 0) {
					new TextDecoder("utf-8").decode(buffer8.slice(old_path_ptr, old_path_ptr + old_path_len));
					new TextDecoder("utf-8").decode(buffer8.slice(new_path_ptr, new_path_ptr + new_path_len));
					return 58;
				} else return 8;
			},
			path_unlink_file(fd, path_ptr, path_len) {
				const buffer8 = new Uint8Array(self.inst.exports.memory.buffer);
				if (self.fds[fd] != void 0) {
					const path = new TextDecoder("utf-8").decode(buffer8.slice(path_ptr, path_ptr + path_len));
					return self.fds[fd].path_unlink_file(path);
				} else return 8;
			},
			poll_oneoff(in_ptr, out_ptr, nsubscriptions) {
				if (nsubscriptions === 0) return 28;
				if (nsubscriptions > 1) {
					debug.log("poll_oneoff: only a single subscription is supported");
					return 58;
				}
				const buffer = new DataView(self.inst.exports.memory.buffer);
				const s = Subscription.read_bytes(buffer, in_ptr);
				const eventtype = s.eventtype;
				const clockid = s.clockid;
				const timeout = s.timeout;
				if (eventtype !== 0) {
					debug.log("poll_oneoff: only clock subscriptions are supported");
					return 58;
				}
				let getNow = void 0;
				if (clockid === 1) getNow = () => BigInt(Math.round(performance.now() * 1e6));
				else if (clockid === 0) getNow = () => BigInt((/* @__PURE__ */ new Date()).getTime()) * 1000000n;
				else return 28;
				const endTime = (s.flags & 1) !== 0 ? timeout : getNow() + timeout;
				while (endTime > getNow());
				new Event(s.userdata, 0, eventtype).write_bytes(buffer, out_ptr);
				return 0;
			},
			proc_exit(exit_code) {
				throw new WASIProcExit(exit_code);
			},
			proc_raise(sig) {
				throw "raised signal " + sig;
			},
			sched_yield() {},
			random_get(buf, buf_len) {
				const buffer8 = new Uint8Array(self.inst.exports.memory.buffer).subarray(buf, buf + buf_len);
				if ("crypto" in globalThis && (typeof SharedArrayBuffer === "undefined" || !(self.inst.exports.memory.buffer instanceof SharedArrayBuffer))) for (let i = 0; i < buf_len; i += 65536) crypto.getRandomValues(buffer8.subarray(i, i + 65536));
				else for (let i = 0; i < buf_len; i++) buffer8[i] = Math.random() * 256 | 0;
			},
			sock_recv(fd, ri_data, ri_flags) {
				throw "sockets not supported";
			},
			sock_send(fd, si_data, si_flags) {
				throw "sockets not supported";
			},
			sock_shutdown(fd, how) {
				throw "sockets not supported";
			},
			sock_accept(fd, flags) {
				throw "sockets not supported";
			}
		};
	}
};
//#endregion
//#region node_modules/@bjorn3/browser_wasi_shim/dist/fd.js
var Inode = class Inode {
	static issue_ino() {
		return Inode.next_ino++;
	}
	static root_ino() {
		return 0n;
	}
	constructor() {
		this.ino = Inode.issue_ino();
	}
};
Inode.next_ino = 1n;
//#endregion
//#region src/wasm-runtime.ts
const encoder = new TextEncoder();
const decoder = new TextDecoder();
var WebPorpidRuntime = class WebPorpidRuntime {
	core;
	constructor(core) {
		this.core = core;
	}
	static async create(module, compiledConfig) {
		const wasi = new WASI([], [], []);
		const instance = await WebAssembly.instantiate(module, { wasi_snapshot_preview1: wasi.wasiImport });
		wasi.initialize(instance);
		const runtime = new WebPorpidRuntime(instance.exports);
		runtime.withOne(compiledConfig, (pointer, length) => runtime.core.wpp_init_config(pointer, length), false);
		return runtime;
	}
	errorText() {
		return decoder.decode(new Uint8Array(this.core.memory.buffer, this.core.wpp_error_ptr(), this.core.wpp_error_len())) || "The webPORPID WASM core failed.";
	}
	resultBytes() {
		return new Uint8Array(this.core.memory.buffer, this.core.wpp_result_ptr(), this.core.wpp_result_len()).slice();
	}
	put(bytes) {
		const pointer = this.core.wpp_alloc(bytes.byteLength);
		if (!pointer && bytes.byteLength) throw new Error("The webPORPID WASM core ran out of memory.");
		new Uint8Array(this.core.memory.buffer, pointer, bytes.byteLength).set(bytes);
		return pointer;
	}
	withOne(bytes, call, result = true) {
		const pointer = this.put(bytes);
		try {
			if (call(pointer, bytes.byteLength) < 0) throw new Error(this.errorText());
			return result ? this.resultBytes() : new Uint8Array();
		} finally {
			this.core.wpp_free(pointer);
		}
	}
	withTwo(first, second, call) {
		const a = this.put(first), b = this.put(second);
		try {
			if (call(a, first.byteLength, b, second.byteLength) < 0) throw new Error(this.errorText());
			return this.resultBytes();
		} finally {
			this.core.wpp_free(a);
			this.core.wpp_free(b);
		}
	}
	preprocess(fastq, firstOrdinal) {
		return this.withOne(encoder.encode(fastq), (pointer, length) => this.core.wpp_preprocess(pointer, length, firstOrdinal));
	}
	partitionCounts(partition) {
		return this.withOne(partition, (pointer, length) => this.core.wpp_partition_counts(pointer, length));
	}
	countFamilies(partition, cutoffs) {
		return this.withTwo(partition, cutoffs, (a, al, b, bl) => this.core.wpp_count_families(a, al, b, bl));
	}
	buildFamilyModel(counts) {
		return this.withOne(counts, (pointer, length) => this.core.wpp_build_family_model(pointer, length));
	}
	initFamilyModel(model) {
		this.withOne(model, (pointer, length) => this.core.wpp_init_family_model(pointer, length), false);
	}
	consensus(partition, cutoffs) {
		return this.withTwo(partition, cutoffs, (a, al, b, bl) => this.core.wpp_consensus_partition(a, al, b, bl));
	}
	stats() {
		if (this.core.wpp_stats() < 0) throw new Error(this.errorText());
		return {
			...JSON.parse(decoder.decode(this.resultBytes())),
			downsampledReads: 0
		};
	}
};
function makeCutoffValues(sampleCounts, maximum) {
	const maxHash = (1n << 64n) - 1n;
	return sampleCounts.map((count) => {
		const cap = maximum < 1 ? count : BigInt(maximum);
		return count === 0n || count <= cap ? maxHash : maxHash * cap / count;
	});
}
function makeCutoffs(sampleCounts, maximum) {
	const writer = new BinaryWriter(), values = makeCutoffValues(sampleCounts, maximum);
	writer.magic("WPT1");
	writer.u32(values.length);
	for (const value of values) writer.u64(value);
	return writer.finish();
}
function decodeFamilyCounts(bytes) {
	const reader = new BinaryReader(bytes);
	reader.magic("WPN1");
	const count = reader.u32();
	const output = [];
	for (let index = 0; index < count; index++) output.push({
		sample: reader.u16(),
		umi: reader.string(),
		count: reader.u32()
	});
	if (!reader.done) throw new Error("Family count payload has trailing bytes.");
	return output;
}
function mergeFamilyCounts(parts) {
	const counts = /* @__PURE__ */ new Map();
	for (const bytes of parts) for (const entry of decodeFamilyCounts(bytes)) {
		const key = `${entry.sample}\0${entry.umi}`, previous = counts.get(key);
		if (previous) previous.count += entry.count;
		else counts.set(key, { ...entry });
	}
	const entries = [...counts.values()].sort((a, b) => a.sample - b.sample || a.umi.localeCompare(b.umi));
	const writer = new BinaryWriter();
	writer.magic("WPN1");
	writer.u32(entries.length);
	for (const entry of entries) {
		writer.u16(entry.sample);
		writer.string(entry.umi);
		writer.u32(entry.count);
	}
	return writer.finish();
}
const DISPOSITIONS = [
	"likely_real",
	"BPB-rejects",
	"heteroduplex",
	"LDA-rejects",
	"UMI_len != 8",
	"family-size-reject"
];
function decodeFamilyModel(bytes, config) {
	const reader = new BinaryReader(bytes);
	reader.magic("WPM1");
	const count = reader.u32();
	const output = [];
	for (let index = 0; index < count; index++) {
		const sampleIndex = reader.u16(), umi = reader.string(), parent = reader.string(), familySize = reader.u32();
		const probability = reader.f64(), disposition = DISPOSITIONS[reader.u8()];
		if (!disposition || !config.samples[sampleIndex]) throw new Error("Family model contains an invalid sample or disposition.");
		output.push({
			sample: config.samples[sampleIndex].name,
			sampleIndex,
			umi: disposition === "BPB-rejects" ? "REJECTED" : umi,
			familySize,
			mostLikelyParent: parent,
			posteriorProbability: probability,
			logOffspringProbability: Math.log(1 - probability),
			disposition
		});
	}
	if (!reader.done) throw new Error("Family model has trailing bytes.");
	return output;
}
function decodeConsensusOutput(bytes, config) {
	const reader = new BinaryReader(bytes);
	reader.magic("WPO1");
	const count = reader.u32();
	const consensuses = [];
	for (let index = 0; index < count; index++) {
		const sampleIndex = reader.u16(), id = reader.string(), umi = reader.string(), familySize = reader.u32();
		const minimumAgreement = reader.f64(), sequence = reader.string(), lowCount = reader.u32();
		const lowAgreementSites = Array.from({ length: lowCount }, () => ({
			position: reader.u32(),
			agreement: reader.f32(),
			modalReadBase: String.fromCharCode(reader.u8()),
			modalRunLength: reader.u32()
		}));
		consensuses.push({
			id,
			sample: config.samples[sampleIndex]?.name ?? String(sampleIndex),
			sampleIndex,
			umi,
			familySize,
			minimumAgreement,
			sequence,
			lowAgreementSites
		});
	}
	const heteroduplexCount = reader.u32(), heteroduplexes = [];
	for (let index = 0; index < heteroduplexCount; index++) heteroduplexes.push(`${reader.u16()}\0${reader.string()}`);
	if (!reader.done) throw new Error("Consensus payload has trailing bytes.");
	return {
		consensuses,
		heteroduplexes
	};
}
function mergeStats(stats, samples) {
	const result = {
		totalReads: 0,
		qualityReads: 0,
		badReads: 0,
		shortReads: 0,
		longReads: 0,
		primerRejects: 0,
		idRejects: 0,
		demultiplexedReads: 0,
		bpbRejects: 0,
		malformedRecords: 0,
		downsampledReads: 0,
		perSample: Array(samples).fill(0)
	};
	for (const part of stats) for (const key of [
		"totalReads",
		"qualityReads",
		"badReads",
		"shortReads",
		"longReads",
		"primerRejects",
		"idRejects",
		"demultiplexedReads",
		"bpbRejects",
		"malformedRecords"
	]) result[key] += part[key];
	for (const part of stats) for (let index = 0; index < part.perSample.length; index++) result.perSample[index] += part.perSample[index];
	return result;
}
//#endregion
export { makeCutoffValues as a, mergeStats as c, decodeFamilyModel as i, WASI as l, decodeConsensusOutput as n, makeCutoffs as o, decodeFamilyCounts as r, mergeFamilyCounts as s, WebPorpidRuntime as t, BinaryWriter as u };
