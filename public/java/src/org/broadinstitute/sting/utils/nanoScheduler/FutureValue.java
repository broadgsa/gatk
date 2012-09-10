package org.broadinstitute.sting.utils.nanoScheduler;

import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;

/**
 * Create a future that simply returns a given value
 *
 * The only standard way to create a future in java is via the ExecutorService interface.
 * If you have a data structure holding futures of value T, and you want to add a
 * value to it for some reason (to add a EOF marker, for instance) you can use this
 * class to create a dummy Future<T> that simply returns a value.
 *
 * @author depristo
 * @since 09/12
 */
class FutureValue<V> implements Future<V> {
    final V value;

    FutureValue(final V value) {
        this.value = value;
    }

    @Override public boolean cancel(boolean mayInterruptIfRunning) {
        return true;
    }

    @Override public boolean isCancelled() {
        return false;
    }

    @Override public boolean isDone() {
        return true;
    }

    @Override public V get() throws InterruptedException, ExecutionException {
        return value;
    }

    @Override public V get(long timeout, TimeUnit unit) throws InterruptedException, ExecutionException, TimeoutException {
        return get();
    }
}
