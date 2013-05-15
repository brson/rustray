extern mod std;

use std::future::{Future,from_port};
use core::comm;

struct ConcurrentCalc<T,U> {
    priv chan: Chan<(T,~fn(T)->U,comm::ChanOne<U>)>
}

impl<T:Owned, U: Owned> ConcurrentCalc<T,U> {
    pub fn new() -> ConcurrentCalc<T,U> {
        let (p, c) = comm::stream();
        do spawn {
            loop {
                match p.try_recv() {
                    Some(message) => {
                        let (data,f,chan): (T,~fn(T)->U,comm::ChanOne<U>) = message;
                        chan.send(f(data))
                    }
                    None => break
                };
            }
        }
        ConcurrentCalc{ chan: c }
    }
    pub fn calculate(&mut self, data: T, f: ~fn(T)->U) -> Future<U> {
        let (p, c) = comm::oneshot();
        self.chan.send( (data,f,c) );
        from_port(p)
    }
}

#[test]
fn testConcurrent() {
    let mut cc: ConcurrentCalc<uint,uint> = ConcurrentCalc::new();
    let mut future = cc.calculate(3, |x| {x+1});
    assert!(future.get() == 4);
    let mut future = cc.calculate(10, |x| {x-4});
    assert!(future.get() == 6);
}

