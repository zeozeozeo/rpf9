//! RPF9 archive reader.
//!
//! The record tables are recovered from the bytes on disk using the cipher and
//! key material embedded from the game executable; the directory is not carved.
//!
//! ```
//! # fn demo(data: &[u8]) -> rpf9::Result<()> {
//! let rpf = rpf9::Rpf::parse(data)?;
//! for entry in rpf.entries() {
//!     let bytes = rpf.read(entry)?;
//!     println!("{:#x} {} bytes", entry.offset, bytes.len());
//! }
//! # Ok(()) }
//! ```
pub mod crypto;
mod entry;
mod keyring;
pub mod reader;

pub mod error;
pub mod header;
#[cfg(feature = "mmap")]
pub mod mmap;

pub const MAGIC: [u8; 4] = *b"9FPR";
pub const HEADER_SIZE: usize = 0x40;

pub const RECORD_SIZE: usize = 16;

pub use entry::Entry;
pub use error::{Error, Result};
pub use header::Header;
#[cfg(feature = "mmap")]
pub use mmap::MappedFile;
pub use reader::Rpf;
